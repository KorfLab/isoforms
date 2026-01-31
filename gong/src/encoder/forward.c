#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"

/* --------------- Forward Algorithm --------------- */
/*
    Yu & Kobayashi 2003, IEEE SPL Vol.10 No.1

    α_t(m, d) = P(O_1...O_t, S_t=m, D_t=d | λ)

    Eq.5:   α_t(m, d) = α_{t-1}(m, d+1) * b_m(O_t)                            [stay]
            α_t(m, d) = [Σ_{n≠m} α_{t-1}(n, 1) * a_{nm}] * b_m(O_t) * p_m(d)  [transition]

    Complexity: O((MD + M²)T)
*/

void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int ekmer = l->B.exon_kmer_len;

    /*
     * INITIALIZATION (Equation 3 in paper):
     *   α_0(m, d) = π_m * b_m(O_1) * p_m(d)
     *
     * For gene prediction, we start in exon state (m=0) with:
     *   - π_exon = 1.0 (always start in exon)
     *   - Accumulate emission probabilities until first potential donor site
     *   - Apply duration probability based on exon length
     */
    if(DEBUG)  printf("Start forward algorithm basis calculation:");

    int     tau;                    // [tau]        :   residential time        aka: possible explicit duration
    int     idx_em_prob, idx_tran_prob;

    double  em_prob;                // [em_prob]    :   bm(o1)                  aka: emission probability
    double  tran_prob;              // [tprob]      :   a(mn)                   aka: transition probability
    double  ed_prob;                // [ed_prob]    :   pm(d)                   aka: explicit duration probability
    double  alpha_sum = LOG_ZERO;   // Initialize to log(0) for log-space accumulation
    // find first donor site
    int  donor = -1;
    char *seq  = info->original_sequence;
    for (int i = FLANK+ed->min_len_exon+1 ; i < info->T-FLANK-ed->min_len_exon-1 ; i++) {
        if (seq[i] == 'G' && seq[i+1] == 'T') {
           donor = i;
           break;
        }
    }

    if (donor == -1) {
        printf("ERROR: No donor site (GT) found in sequence!\n");
        printf("Search range: %d to %d\n", FLANK+ed->min_len_exon, info->T-FLANK-ed->min_len_exon);
        printf("\n");
        return;
    }

    alpha->first_dons = donor;

    // boundary check
    tau = info->T-FLANK-ed->min_len_exon - donor;
    if (tau > ed->max_len_exon)   tau = ed->max_len_exon;

    if (DEBUG)  printf("Processing exon region from %d to %d (tau=%d)\n", FLANK, donor, tau);

    // get emission prob before first donor site
    alpha->basis[0][0] = LOG_ZERO;
    alpha_sum = LOG_ZERO;  // Initialize to log(0)

    int duration = -1;
    for (int bps = FLANK ; bps < donor+1 ; bps++) {
        // Boundary check: ensure we don't go negative
        int kmer_start = bps - ekmer + 1;
        if (kmer_start < 0) continue;

        idx_em_prob = base4_to_int(info->numerical_sequence, kmer_start, ekmer);
        em_prob     = l->B.exon[idx_em_prob];

        // Log-space accumulation
        if (IS_LOG_ZERO(alpha_sum)) alpha_sum = em_prob;
        else if (!IS_LOG_ZERO(em_prob)) alpha_sum = alpha_sum + em_prob;

        duration++;

        if (bps == donor-1 || bps == donor) {
            ed_prob = ed->exon[duration];
            alpha->a[bps][0] = LOG_MUL(alpha_sum, ed_prob);
        }

        if (bps == donor) {
            alpha->basis[0][0] = alpha_sum;

            for (int i = duration ; i < tau - duration ; i++) {
                ed_prob = ed->exon[i];
                alpha->basis[0][i-duration] = LOG_MUL(alpha_sum, ed_prob);
            }
        }
    }

    if (DEBUG)  printf("Exon basis calculation complete. alpha_sum = %e\n", alpha_sum);
    /*
        for intron basis
        we need wait until the first donor site appear for continue calculation
            a(0)(intron, d) = α(-1)(exon, 0) * a(exon|intron) * bintron(o0) * pintron(d)
    aka               total = pi * tprob * emprob * edprob
    */
    if  (tau > ed->max_len_intron)  tau = ed->max_len_intron;

    em_prob         = l->B.intron[idx_em_prob];
    idx_tran_prob   = base4_to_int(info->numerical_sequence, donor, dkmer);
    tran_prob       = l->A.dons[idx_tran_prob];

    // Log-space: multiply is addition, check for LOG_ZERO
    if (IS_LOG_ZERO(tran_prob) || IS_LOG_ZERO(alpha->a[donor-1][0])) {
        alpha_sum = LOG_ZERO;
    } else {
        alpha_sum = tran_prob + alpha->a[donor-1][0] + em_prob;
    }

    for (int d = 0 ; d < tau ; d++) {
        ed_prob = ed->intron[d];
        alpha->basis[1][d] = LOG_MUL(alpha_sum, ed_prob);
    }

    alpha->a[donor][1] = alpha->basis[1][0];

    if(DEBUG) {
        printf("=== Basis for forward Algorithm Debug ===\n");
        printf("Donor: %1e, %1e", alpha->a[donor-1][0], alpha->a[donor][0]);
        printf("Acceptor:%1e, %1e", alpha->a[donor-1][1], alpha->a[donor][1]);
    }

    if (DEBUG)      printf("\tFinished\n");
}

/*
    Forward Recursion
        basis[m][d] : accumulated forward prob for state m, remaining duration d
        a[t][m]     : forward prob at position t, state m, duration=1

    Transition via splice signal PWM:
        exon→intron : donor PWM at GT
        intron→exon : acceptor PWM at AG
*/
void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    int ikmer = l->B.intron_kmer_len;
    if(DEBUG)  printf("Start computation for forward algorithm:\n");

    int start = alpha->first_dons+1;
    int end   = info->T-FLANK-ed->min_len_exon;
    int tau   = end-start+1;
    int bound;

    for (int bps = start; bps < end; bps++) {
        tau--;

        for (int hs = 0; hs < HS; hs++) {
            double  tran_prob, ed_prob, em_prob, tran_node, cont_node;
            int     idx_tran_prob, idx_em_prob, con_hs, okmer;
            // observed sequence kmer
            okmer = (hs == 0) ? ekmer : ikmer;
            bound = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;
            if (tau >= bound)   tau = bound;

            con_hs      = (hs == 0) ? 1 : 0;
            idx_em_prob = base4_to_int(info->numerical_sequence, bps-okmer+1, okmer);
            em_prob     = (hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];

            if (hs == 0) {
                idx_tran_prob = base4_to_int(info->numerical_sequence, bps-akmer, akmer);
                tran_prob = l->A.accs[idx_tran_prob];
            } else {
                idx_tran_prob = base4_to_int(info->numerical_sequence, bps, dkmer);
                tran_prob = l->A.dons[idx_tran_prob];
            }
            tran_node = alpha->a[bps - 1][con_hs];

            for (int i = 0; i < tau; i++) {
                cont_node = (i != tau - 1) ? alpha->basis[hs][i+1] : LOG_ZERO;

                double cont_sum, tran_sum;

                // Continuation: stay in same state, duration decreases
                if (IS_LOG_ZERO(cont_node) || IS_LOG_ZERO(em_prob)) {
                    cont_sum = LOG_ZERO;
                } else {
                    cont_sum = cont_node + em_prob;
                }

                ed_prob = (hs == 0) ? ed->exon[i] : ed->intron[i];

                // Transition: come from other state
                if (IS_LOG_ZERO(tran_node) || IS_LOG_ZERO(tran_prob) || IS_LOG_ZERO(ed_prob)) {
                    tran_sum = LOG_ZERO;
                } else {
                    tran_sum = tran_node + tran_prob + ed_prob + em_prob;
                }

                // Combine continuation and transition using log-sum-exp
                if (IS_LOG_ZERO(tran_sum)) {
                    alpha->basis[hs][i] = cont_sum;
                } else if (IS_LOG_ZERO(cont_sum)) {
                    alpha->basis[hs][i] = tran_sum;
                } else {
                    alpha->basis[hs][i] = log_add(cont_sum, tran_sum);
                }
                alpha->a[bps][hs] = alpha->basis[hs][0];
            }
        }
    }

    if(DEBUG)  printf("\tFinished\n");
}
