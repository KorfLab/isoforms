#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"

/* --------------- Backward Algorithm --------------- */
/*
    Yu & Kobayashi 2003, Eq.8

    β_t(m, d) = P(O_{t+1}...O_T | S_t=m, D_t=d, λ)

    Eq.8:   β_t(m, d) = b_m(O_{t+1}) * β_{t+1}(m, d-1)                           [stay]
            β_t(m, 1) = Σ_{n≠m} a_{mn} * b_n(O_{t+1}) * Σ_{d'≥1} p_n(d') * β_{t+1}(n, d')  [transition]

    Initialize: β_T(m, d) = 1 for all m, d
*/

void basis_bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    if(DEBUG)  printf("Start Backward Algorithm basis calculation:\n");

    double em_prob;             // [em_prob]    : emission probability
    double ed_prob;             // [ed_prob]    : explicit duration probability
    double bw_sum = LOG_ZERO;   // [bw_sum]     : arbitrary reused computation variable (log-space)
    double tran_prob;           // [tran_prob]  : a(mn)

    int    idx_em_prob, idx_tran_prob;
    // find first accs site from RHS
    char *seq = info->original_sequence;
    int  accs = -1;  // Initialize to invalid value
    for (int i = info->T-FLANK-ed->min_len_exon-1 -1; i > FLANK+ed->min_len_exon; i--) {
        if (seq[i-1] == 'A' && seq[i] == 'G') {
           accs = i+1;
           break;
       }
    }

    // Error check: no acceptor site found
    if (accs == -1) {
        printf("ERROR: No acceptor site (AG) found in sequence!\n");
        printf("Search range: %d to %d\n", FLANK+ed->min_len_exon, info->T-FLANK-ed->min_len_exon);
        beta->last_accs = FLANK + ed->min_len_exon;  // Set to safe default
        return;
    }

    if (DEBUG)  printf("\nStart compute backward basis from %d to %d", accs, info->T-FLANK);
    beta->last_accs = accs;
    // compute the bw exon until first accs site
    bw_sum = LOG_ZERO;
    for (int bps = info->T-FLANK ; bps > accs; bps--) {
        // Boundary check
        int kmer_start = bps - ekmer + 1;
        if (kmer_start < 0) continue;

        idx_em_prob = base4_to_int(info->numerical_sequence, kmer_start, ekmer);
        em_prob     = l->B.exon[idx_em_prob];

        // Log-space accumulation
        if (IS_LOG_ZERO(bw_sum)) bw_sum = em_prob;
        else if (!IS_LOG_ZERO(em_prob)) bw_sum = bw_sum + em_prob;
    }

    int exon_len = info->T-FLANK-accs;
    if (exon_len >= 0 && exon_len < ed->exon_len) {
        beta->basis[0][exon_len] = bw_sum;
    }

    // update the intron basis
    idx_tran_prob  = base4_to_int(info->numerical_sequence, accs-akmer+1, akmer);
    tran_prob      = l->A.accs[idx_tran_prob];
    ed_prob        = (exon_len >= 0 && exon_len < ed->exon_len) ? ed->exon[exon_len] : LOG_ZERO;

    if (IS_LOG_ZERO(tran_prob) || IS_LOG_ZERO(bw_sum)) {
        bw_sum = LOG_ZERO;
    } else {
        bw_sum = bw_sum + ed_prob + tran_prob;
    }

    // for intron | exon(min_exon_len) ; FLANK
    beta->basis[1][0] = bw_sum;
    beta->b[accs][1]  = bw_sum;

    if(DEBUG)  printf("\tFinished\n");
}

/*
    Backward Recursion
        basis[m][d] : accumulated backward prob for state m, remaining duration d
        b[t][m]     : backward prob at position t, state m

    Process from T→1, symmetric to forward pass
*/
void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    int ikmer = l->B.intron_kmer_len;
    if (DEBUG)      printf("Start Backward Algorithm:");

    int start = beta->last_accs-1;
    int end   = FLANK;
    int tau;                    // [tau]    : the residual bound for hidden state

    if (DEBUG)    printf("\n\tCompute from region %d to %d", start, end);

    for (int bps = start ; bps > end ; bps--) {
        double ed_prob;         // [ed_prob]    : pn(d)                                 aka: explicit duration probability
        double node;
        double tran_prob;       // [tran_prob]  : a(mn)                                 aka: transition probability
        double em_prob;         // [em_prob]    : bn(ot+1)                              aka: emission probability
        double bw_sum;          // [bw_sum]     : see note in first part                aka: backward sum

        int    con_hs;          // [con_hs]     : conjugated hidden state
        int    idx_tran_prob, idx_em_prob, okmer;

        for (int hs = 0 ; hs < HS ; hs++) {

            con_hs = (hs == 0) ? 1 : 0 ;
            tau    = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;
            // make summation for all previous layer of node (log-sum-exp over durations)
            int n = 0;
            for (int i = 0 ; i < tau ; i++) {
                node = beta->basis[con_hs][i];
                if (IS_LOG_ZERO(node)) continue;

                ed_prob = (con_hs == 0) ? ed->exon[i] : ed->intron[i];
                if (IS_LOG_ZERO(ed_prob)) continue;

                l->log_values[n] = node + ed_prob;
                n++;
            }
            bw_sum = log_sum_exp(l->log_values, n);

            if (hs == 0) {
                idx_tran_prob = base4_to_int(info->numerical_sequence, bps+1, dkmer);
                tran_prob     = l->A.dons[idx_tran_prob];
            } else {
                idx_tran_prob = base4_to_int(info->numerical_sequence, bps-akmer+1, akmer);
                tran_prob     = l->A.accs[idx_tran_prob];
            }

            okmer       = (con_hs == 0) ? ekmer : ikmer;
            // Note: Using bps-okmer+1 to match forward algorithm indexing
            idx_em_prob = base4_to_int(info->numerical_sequence, bps-okmer+1, okmer);
            em_prob     = (con_hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];

            if (IS_LOG_ZERO(tran_prob) || IS_LOG_ZERO(bw_sum)) {
                bw_sum = LOG_ZERO;
            } else {
                bw_sum = tran_prob + bw_sum + em_prob;
            }

            beta->b[bps][hs] = bw_sum;
        }

        for (int hs = 0 ; hs < HS ; hs++) {
            okmer       = (hs == 0) ? ekmer : ikmer;
            idx_em_prob = base4_to_int(info->numerical_sequence, bps-okmer+1, okmer);
            em_prob     = (hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];
            tau         = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;
            for (int i = tau-1 ; i > 0 ; --i) {
                node = beta->basis[hs][i-1];
                if (IS_LOG_ZERO(node) || IS_LOG_ZERO(em_prob)) {
                    beta->basis[hs][i] = LOG_ZERO;
                } else {
                    beta->basis[hs][i] = node + em_prob;
                }
            }
            // update the b[t][0]
            beta->basis[hs][0] = beta->b[bps][hs];
        }
    }
    if(DEBUG)  printf("\tFinished.\n");
}
