#include <stdio.h>
#include <stdlib.h>
#include "model.h"

/* --------------- Posterior Probability --------------- */

void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Pos_prob *pos) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    /*
        Posterior probability (Yu-Kobayashi paper):
            ξ_t(m, n) = α_{t-1}(m, 1) * a_{mn} * b_n(O_t) * Σ_{d≥1} p_n(d) * β_t(n, d)

        Simplified (since we store α(t)(m,1) and β(t)(m,1) with transitions included):
            ξ(t)(m, n) ≈ α(t-1)(m, 1) * β(t-1)(m, 1)

        In log-space: xi = fw + bw (both are already in log-space)
    */

    double fw;      // [fw]    :   α(t-1)(m, 1) in log-space
    double bw;      // [bw]    :   β(t)(m, 1) in log-space
    double xi;      // [xi]    :   posterior probability in log-space

    for (int bps = FLANK ; bps < info->T-FLANK ; bps++) {
        // Donor site posterior: P(transition exon→intron at position bps)
        fw = alpha->a[bps-1][0];  // Forward prob of being in exon at bps-1
        bw = beta->b[bps-1][0];   // Backward prob from exon at bps-1

        if (!IS_LOG_ZERO(fw) && !IS_LOG_ZERO(bw)) {
            xi = fw + bw;  // Log-space multiplication
            pos->xi[bps][0] = xi;
        } else {
            pos->xi[bps][0] = LOG_ZERO;
        }

        // Acceptor site posterior: P(transition intron→exon at position bps)
        fw = alpha->a[bps][1];    // Forward prob of being in intron at bps
        bw = beta->b[bps][1];     // Backward prob from intron at bps

        if (!IS_LOG_ZERO(fw) && !IS_LOG_ZERO(bw)) {
            xi = fw + bw;  // Log-space multiplication
            pos->xi[bps][1] = xi;
        } else {
            pos->xi[bps][1] = LOG_ZERO;
        }
    }
}

void parse_splice_sites(Pos_prob *pos, Observed_events *info) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;

    // Count donor sites (non-LOG_ZERO posterior probabilities)
    int dons = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (!IS_LOG_ZERO(pos->xi[i][0])) dons++;
    }
    pos->dons_val = malloc(dons * sizeof(double));
    pos->dons_bps = malloc(dons * sizeof(int));
    pos->dons     = dons;

    int idx = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (!IS_LOG_ZERO(pos->xi[i][0])) {
            pos->dons_val[idx] = pos->xi[i][0];
            pos->dons_bps[idx] = i;
            idx++;
        }
    }

    // Count acceptor sites
    int accs = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (!IS_LOG_ZERO(pos->xi[i][1])) accs++;
    }
    pos->accs_val = malloc(accs * sizeof(double));
    pos->accs_bps = malloc(accs * sizeof(int));
    pos->accs     = accs;
    idx = 0;
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (!IS_LOG_ZERO(pos->xi[i][1])) {
            pos->accs_val[idx] = pos->xi[i][1];
            pos->accs_bps[idx] = i;
            idx++;
        }
    }
}

void free_splice_sites(Pos_prob *pos) {
    if (pos->dons_val) free(pos->dons_val);
    if (pos->dons_bps) free(pos->dons_bps);
    if (pos->accs_val) free(pos->accs_val);
    if (pos->accs_bps) free(pos->accs_bps);
}
