#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "model.h"

/* --------------- Initializer --------------- */

void numerical_transcription(Observed_events *info, const char *seq) {
    if(DEBUG)  printf("Start transforming original sequence into base4:\n");
    // turns original sequence into int 
    size_t len  = strlen(seq);
    info->T     = len;
    info->numerical_sequence = malloc ( len * sizeof(int) );

    for (size_t i = 0; i < len; i++) {
        if      (seq[i] == 'A')     info->numerical_sequence[i] = 0;
        else if (seq[i] == 'C')     info->numerical_sequence[i] = 1;
        else if (seq[i] == 'G')     info->numerical_sequence[i] = 2;
        else if (seq[i] == 'T')     info->numerical_sequence[i] = 3;
    }

    if(DEBUG)  printf("\tWe get numerical sequence with Seq len: %d\n", info->T);
    if(DEBUG)  printf("\tFinished\n");
    if(DEBUG)  printf("\n");
}

/* --------------- Computation Function --------------- */

int power(int base, int exp) {                                          // wtf, C don't have power for int
    int result = 1;
    for (int i = 0 ; i < exp ; i++) {
        result *= base;
    }
    return result;
}

int base4_to_int(int *array, int beg, int length) {
    int values = 0;
    int value;
    for (int i = 0; i < length; i++) {
        value  =  array[beg + i];
        values += value * power(4, length - i - 1);
    }
    return values;
}

double log_sum_exp(double *array, int n) {
    if(!array){
        printf("Something wrong with log_sum_exp trick. Invalid Input\n");
        return 0.0;
    }
    if(n <= 0) return 0.0;

    double max = array[0];
    double sum = 0.0;
    // find max
    for (int i = 1 ; i < n ; i++) {
        if( array[i] > max )     max = array[i];                  
    }
    // computation
    for (int i = 0 ; i < n ; i++) {
        sum += exp(array[i] - max);
    }
    return max + log(sum);       
}

/* --------------- Validator --------------- */

void tolerance_checker(double *array, int len, const double epsilon) {
    if(!array){
        printf("This is not a valid input for float20.0\n");
        return;
    }

    for(int i = 0 ; i < len ; i++) {
        if( array[i] < epsilon )    array[i] = 0;
    }
}

void log_space_converter(double *array, int len) {
    if(!array){
        printf("This is not a valid input for log_space converter\n");
        return;
    }

    for (int i = 0; i < len ; i++) {
        if( array[i] == 0.0 )   continue;
        else                    array[i] = log(array[i]);
    }
}

/* --------------- Forward Algorithm --------------- */

void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    int ikmer = l->B.intron_kmer_len;
    /*
        given initial formula
            a(0)(m, d) = pi(m) * bm(d) * pm(d)

        first part for exon basis
        notice:         since initial min_len_exon bound are 100% exon; we need get there emission prob product
                        correct explicit duration probability
    in our case         = ∏(0 -> min_len_exon)bm(d) * pm(d+min_len_exon)
    aka:          total = pi * ed_prob
    */
    if(DEBUG)  printf("Start forward algorithm basis calculation:");

    int     tau;                    // [tau]        :   residential time        aka: possible explicit duration
    int     idx_em_prob, idx_tran_prob;

    double  em_prob;                // [em_prob]    :   bm(o1)                  aka: emission probability
    double  tran_prob;              // [tprob]      :   a(mn)                   aka: transition probability
    double  ed_prob;                // [ed_prob]    :   pm(d)                   aka: explicit duration probability
    double  alpha_sum = 0.0;
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
    alpha->basis[0][0] = 0.0;

    int duration = -1;
    for (int bps = FLANK ; bps < donor+1 ; bps++) {
        idx_em_prob = base4_to_int(info->numerical_sequence, bps-ekmer+1, ekmer);
        em_prob     = l->B.exon[idx_em_prob];
        alpha_sum   = alpha_sum+em_prob;
        duration++;

        if (bps == donor-1 || bps == donor) {
            ed_prob      = ed->exon[duration];

            if (ed_prob == 0.0) alpha->a[bps][0]    = 0.0;
            else                alpha->a[bps][0]    = alpha_sum+ed_prob;
        }

        if (bps == donor) {   
            alpha->basis[0][0] = alpha_sum;

            for (int i = duration ; i < tau - duration ; i++) {
                ed_prob = ed->exon[i];
                if (ed_prob == 0.0) alpha->basis[0][i-duration] = 0.0;
                else                alpha->basis[0][i-duration] = alpha_sum+ed_prob;
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

    if (tran_prob == 0.0)   alpha_sum = 0.0;
    else                    alpha_sum = tran_prob+alpha->a[donor-1][0]+em_prob;

    for (int d = 0 ; d < tau ; d++) {
        ed_prob = ed->intron[d];

        if  (ed_prob == 0.0)    alpha->basis[1][d] = 0.0;
        else                    alpha->basis[1][d] = alpha_sum+ed_prob;
    }

    alpha->a[donor][1] = alpha->basis[1][0]; 

    if(DEBUG) {
        printf("=== Basis for forward Algorithm Debug ===\n");
        printf("Donor: %1e, %1e", alpha->a[donor-1][0], alpha->a[donor][0]);
        printf("Acceptor:%1e, %1e", alpha->a[donor-1][1], alpha->a[donor][1]);
    }

    if (DEBUG)      printf("\tFinished\n");
}

void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    int ikmer = l->B.intron_kmer_len;
    /*
        Forward recursion formula for duration-based HMM:

        α_t(m, d) =   α_(t-1)(m, d+1) * b_m(o_t)
                + [sum over n ≠ m of α_(t-1)(n, 1) * a_nm] * b_m(o_t) * p_m(d)

        Logic:
            1. If you're still in the same state (m), just extend duration: 
                   take α from previous time, same state, one longer duration.
            2. Or you just transitioned into state m from another state n:
                   sum over all other states n, and get their α with duration 1,
                   multiply by transition probability from n to m,
                   then add emission and duration prob for new state m.
    */
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
                cont_node = (i != tau - 1) ? alpha->basis[hs][i+1] : 0.0;

                double cont_sum, tran_sum;

                if (cont_node == 0.0)   cont_sum = 0.0;
                else                    cont_sum = cont_node + em_prob;

                ed_prob = (hs == 0) ? ed->exon[i] : ed->intron[i];

                if (tran_node == 0.0 || tran_prob == 0.0 || ed_prob == 0.0) {
                    tran_sum = 0.0;
                } else {
                    tran_sum = tran_node+tran_prob+ed_prob+em_prob;
                }

                if      (tran_sum == 0.0)   alpha->basis[hs][i] = cont_sum;
                else if (cont_sum == 0.0)   alpha->basis[hs][i] = tran_sum;
                else {
                    double logs[2] = { cont_sum, tran_sum };
                    double total   = log_sum_exp(logs, 2);
                    alpha->basis[hs][i] = total;
                }
                alpha->a[bps][hs] = alpha->basis[hs][0];
            }
        }
    }

    if(DEBUG)  printf("\tFinished\n");
}

/* --------------- Backward Algorithm --------------- */

void basis_bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    int ikmer = l->B.intron_kmer_len;
    /*   
        for the most right bound
        given initial condition
            β T(m, d) = 1
            β t(m, d) = bm(ot+1)*β T(m, d-1)
    in short       pi = ∏(min_len_exon - 1) emission probability (lot of logic skip here, refer paper)

        computation
            β t(m, 1) = amn * bn(Ot + 1) * Σ(d>=1) pn(d) * β(t + 1)(n , d)
    aka         total = tprob * pi

        boundary
            calculate until we reach first donor site
    */
    if(DEBUG)  printf("Start Backward Algorithm basis calculation:\n");

    double em_prob;             // [em_prob]    : emission probability
    double ed_prob;             // [ed_prob]    : explicit duration probability
    double bw_sum = 0.0;        // [bw_sum]     : arbitrary reused computation variable
    double tran_prob;           // [tran_prob]  : a(mn)

    int    idx_em_prob, idx_tran_prob;
    // find first accs site from RHS
    char *seq = info->original_sequence;
    int  accs;
    for (int i = info->T-FLANK-ed->min_len_exon-1 -1; i > FLANK+ed->min_len_exon; i--) {
        if (seq[i-1] == 'A' && seq[i] == 'G') {
           accs = i+1;
           break;
       }
    }

    if (DEBUG)  printf("\nStart compute backward basis from %d to %d", accs, info->T-FLANK);
    beta->last_accs = accs;
    // compute the bw exon until first accs site
    for (int bps = info->T-FLANK ; bps > accs; bps--) {
        idx_em_prob = base4_to_int(info->numerical_sequence, bps-ekmer+1, ekmer);
        em_prob     = l->B.exon[idx_em_prob];
        bw_sum      = bw_sum+em_prob;
    }

    int exon_len = info->T-FLANK-accs;
    beta->basis[0][exon_len] = bw_sum;
    // update the intron basis
    idx_tran_prob  = base4_to_int(info->numerical_sequence, accs-akmer+1, akmer);
    tran_prob      = l->A.accs[idx_tran_prob];
    ed_prob        = ed->exon[exon_len];

    if (tran_prob == 0.0)   bw_sum = 0.0;
    else                    bw_sum = bw_sum+ed_prob+tran_prob;

    // for intron | exon(min_exon_len) ; FLANK
    beta->basis[1][0] = bw_sum;
    beta->b[accs][1]  = bw_sum;

    if(DEBUG)  printf("\tFinished\n");
}

void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dkmer = l->B.don_kmer_len;
    int akmer = l->B.acc_kmer_len;
    int ekmer = l->B.exon_kmer_len;
    int ikmer = l->B.intron_kmer_len;
    /*
    Backward recursion formula for duration-based HMM:

    Case 1: d > 1
    β_t(m, d) = b_m(o_{t+1}) * β_{t+1}(m, d - 1)

        Logic:
        - When remaining duration d > 1, we stay in the same state m.
        - Emit observation o_{t+1} using b_m(o_{t+1}).
        - Reduce duration by 1 and move to next time step.

    Case 2: d == 1
    β_t(m, 1) = sum over n ≠ m of [ a_mn * b_n(o_{t+1}) * sum_{d' ≥ 1} (p_n(d') * β_{t+1}(n, d')) ]

        Logic:
        - When duration is exhausted (d == 1), we must transition to a new state n ≠ m.
        - Transition with a_mn, emit o_{t+1} with b_n(o_{t+1}),
          and sum over all durations d' ≥ 1 of the new state with duration probability p_n(d').
    */
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
            // make summation for all previous layer of node
            int n = 0;
            for (int i = 0 ; i < tau ; i++) {
                node = beta->basis[con_hs][i];
                if (node == 0.0)    continue;
                
                ed_prob = (con_hs == 0) ? ed->exon[i] : ed->intron[i];
                if (ed_prob == 0.0)     continue;
                n++;
                l->log_values[n-1] = node+ed_prob;
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
            idx_em_prob = base4_to_int(info->numerical_sequence, bps-okmer+2, okmer);
            em_prob     = (con_hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];

            if (tran_prob == 0.0 || bw_sum == 0.0)  bw_sum = 0.0;
            else                                    bw_sum = tran_prob+bw_sum+em_prob;

            beta->b[bps][hs] = bw_sum;
        }

        for (int hs = 0 ; hs < HS ; hs++) {
            okmer       = (hs == 0) ? ekmer : ikmer;
            idx_em_prob = base4_to_int(info->numerical_sequence, bps-okmer+2, okmer);
            em_prob     = (hs == 0) ? l->B.exon[idx_em_prob] : l->B.intron[idx_em_prob];
            tau         = (hs == 0) ? ed->max_len_exon : ed->max_len_intron;
            for (int i = tau-1 ; i > 0 ; --i) {
                node = beta->basis[hs][i-1];
                beta->basis[hs][i] = (node == 0.0) ? 0.0 : node+em_prob;
            } 
            // update the b[t][0]
            beta->basis[hs][0] = beta->b[bps][hs];
        }
    }
    if(DEBUG)  printf("\tFinished.\n");
}

/* --------------- Posterior Probability --------------- */

void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Pos_prob *pos) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    /*
        update the posterior probability coupled inside bw algo
            ξ(t)(m, n) = α(t-1)(m, 1) * a(mn) * bn(ot) * Σ(n)bm(ot) * ed(prob)
        aka
            ξ(t)(m, n) = α(t-1)(m, 1) * β(t)(m, 1)
        which is
            xi = fw * bw
    */
    
    double fw;      // [fw]    :   α(t-1)(m, 1)
    double bw;      // [bw]    :   β(t)(m, 1)
    double xi;      // [xi]    :   greek letter for pos_prob
    
    for (int bps = FLANK ; bps < info->T-FLANK ; bps++) {
        //donor
        fw = alpha->a[bps-1][0];
        bw = beta->b[bps-1][0];
        if (fw != 0.0 && bw != 0.0) {
            xi = fw+bw;
            pos->xi[bps][0] = xi;
        } else {
            pos->xi[bps][0] = 0.0;
        }

        // accetor
        fw = alpha->a[bps][1];
        bw = beta->b[bps][1];
        if (fw != 0.0 && bw != 0.0) {
            xi = fw+bw;
            pos->xi[bps][1] = xi;
        } else {
            pos->xi[bps][1] = 0.0;
        }
    }
}

void parse_splice_sites(Pos_prob *pos, Observed_events *info) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int dons = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][0] != 0.0) dons++;
    }
    pos->dons_val = malloc(dons * sizeof(double));
    pos->dons_bps = malloc(dons * sizeof(int));
    pos->dons     = dons;

    int idx = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][0] != 0.0) {
            pos->dons_val[idx] = pos->xi[i][0];
            pos->dons_bps[idx] = i;
            idx++;
        }
    }

    int accs = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][1] != 0.0) accs++;
    }
    pos->accs_val = malloc(accs * sizeof(double));
    pos->accs_bps = malloc(accs * sizeof(int));
    pos->accs     = accs;
    idx = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][1] != 0.0) {
            pos->accs_val[idx] = pos->xi[i][1];
            pos->accs_bps[idx] = i;
            idx++;
        }
    }
}

/* --------------- Memory Allocation --------------- */

void allocate_fw(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed) {
    if(DEBUG)  printf("Start allocate memory for the forward algorithm:");

    // alpha->a[t][i]
    int size_array = info->T;
    alpha->a = malloc ( (size_array) * sizeof(double*) );           // [t]: specific time among all observed events 

    for( int i = 0 ; i < size_array; i++ )                          // [0] for exon ; [1] for intron
        alpha->a[i] = calloc( HS , sizeof(double) );                // each spot is storing a(t)(m, 1)           
    // alpha->basis[i][d]
    alpha->basis    = malloc( HS * sizeof(double*) );               // [0] for exon ; [1] for intron
    // max duration for exon or intron
    alpha->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    alpha->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    if(DEBUG)  printf("\tFinished\n");
}

void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info) {
    if(DEBUG)  printf("Start allocate memory for the backward algorithm:");
                                    
    // allocate basis
    beta->basis    = malloc( HS * sizeof(double*) );                    // [0] for exon ; [1] for intron
    beta->basis[0] = calloc( ed->max_len_exon  , sizeof(double) );      // max duration for exon or intron
    beta->basis[1] = calloc( ed->max_len_intron, sizeof(double) );

    // allocate storage array
    int size_array = info->T;
    beta->b = malloc( (size_array) * sizeof(double*) );                 // [0] for exon ; [1] for intron
    for( int i = 0 ; i < size_array; i++ )                             
        beta->b[i] = calloc( HS , sizeof(double) );  

    if(DEBUG)  printf("\tFinished\n");
}

void allocate_pos(Pos_prob *pos, Observed_events *info){
    if(DEBUG)  printf("Start Initialize Data Structure for Posterior Probability");

    int sarray  = info->T;
    pos->xi     = malloc (sarray * sizeof(double*)); 
    for (int i = 0 ; i < sarray; i++)
        pos->xi[i] = calloc( HS , sizeof(double) );       
    
    if(DEBUG)  printf("\tFinished\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha) {
    if(DEBUG)  printf("Clearing up forward algorithm memory:\n");

    int size_array = info->T;
    for (int i = 0; i < size_array; i++) 
        free(alpha->a[i]);

    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    if(DEBUG)  printf("\tFinished\n");
}

void free_beta(Observed_events *info, Backward_algorithm *beta) {
    if(DEBUG)  printf("Clearing up backward algorithm memory:\n");

    int size_array = info->T;
    for (int i = 0; i < size_array; i++)
        free(beta->b[i]);

    free(beta->b);
    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    if(DEBUG)  printf("\tFinished\n");
}

void free_pos(Pos_prob *pos, Observed_events *info) {
    if(DEBUG)  printf("Clearing up Posterior Probabilities:\n");

    int T = info->T;
    for (int t = 0; t < T; t++)
        free(pos->xi[t]);
    free(pos->xi);

    if(DEBUG)  printf("\tFinished\n");
}