#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"
#include "randomf.h"

/* --------------- Auxilary Function --------------- */

static double log_sum_sub(double val, double add, double sub) {
    double max          = (val > add) ? val : add;    
    double sum          = exp(val - max) + exp(add - max);
    double log_sum      = max + log(sum);
    double max_val      = (log_sum > sub) ? log_sum : sub;    
    double diff         = exp(log_sum - max_val) - exp(sub - max_val);

    if (diff <= 0) {
        return NAN;
    }

    return max_val + log(diff);
}

/* --------------- Local Viterbi Allocation --------------- */

static Vitbi_algo* allocate_local_viterbi(Observed_events *info) {
    Vitbi_algo *vit = malloc(sizeof(Vitbi_algo));
    vit->v = malloc(HS * sizeof(double*));
    for(int i = 0; i < HS; i++) {
        vit->v[i] = calloc(info->T, sizeof(double));
    }
    return vit;
}

static void free_local_viterbi(Vitbi_algo *vit) {
    for(int i = 0; i < HS; i++) {
        free(vit->v[i]);
    }
    free(vit->v);
    free(vit);
}

/* --------------- Validation Function --------------- */

int validate_isoform(Isoform *iso, Explicit_duration *ed) {
    // Check if there are introns
    if (iso->n_introns == 0) {
        return 1;  // No introns is valid (single exon gene)
    }
    
    // Check that we have matching donors and acceptors
    if (iso->dons == NULL || iso->accs == NULL) {
        return 0;
    }
    
    // Validate each intron
    for (int i = 0; i < iso->n_introns; i++) {
        // Donor must come before acceptor
        if (iso->dons[i] >= iso->accs[i]) {
            return 0;  // Invalid: donor after or at acceptor
        }
        
        // Check minimum intron length
        int intron_length = iso->accs[i] - iso->dons[i];
        if (intron_length < ed->min_len_intron) {
            return 0;  // Invalid: intron too short
        }
        
        // Check maximum intron length
        if (intron_length > ed->max_len_intron) {
            return 0;  // Invalid: intron too long
        }
        
        // If not the first intron, check exon between previous acceptor and current donor
        if (i > 0) {
            int exon_length = iso->dons[i] - iso->accs[i-1];
            if (exon_length < ed->min_len_exon) {
                return 0;  // Invalid: exon too short
            }
            if (exon_length > ed->max_len_exon) {
                return 0;  // Invalid: exon too long
            }
        }
    }
    
    // Check first exon (before first donor)
    if (iso->n_introns > 0) {
        int first_exon_length = iso->dons[0] - iso->beg;
        if (first_exon_length < ed->min_len_exon) {
            return 0;  // Invalid: first exon too short
        }
    }
    
    // Check last exon (after last acceptor)
    if (iso->n_introns > 0) {
        int last_exon_length = iso->end - iso->accs[iso->n_introns - 1];
        if (last_exon_length < ed->min_len_exon) {
            return 0;  // Invalid: last exon too short
        }
    }
    
    return 1;  // Valid isoform
}

/* --------------- Debug Function --------------- */

void debug_isoform(Isoform *iso, const char *context) {
    printf("\n=== Isoform Debug (%s) ===\n", context);
    printf("Begin: %d, End: %d\n", iso->beg, iso->end);
    printf("Number of introns: %d\n", iso->n_introns);
    
    if (iso->n_introns > 0) {
        printf("Splice sites:\n");
        for (int i = 0; i < iso->n_introns; i++) {
            int intron_len = iso->accs[i] - iso->dons[i];
            printf("  Intron %d: Donor=%d, Acceptor=%d, Length=%d", 
                   i+1, iso->dons[i], iso->accs[i], intron_len);
            
            // Flag problematic introns
            if (iso->dons[i] >= iso->accs[i]) {
                printf(" [ERROR: Donor >= Acceptor!]");
            } else if (intron_len < 20) {  // Assuming min intron is ~20bp
                printf(" [WARNING: Very short intron!]");
            }
            printf("\n");
            
            // Check exon between introns
            if (i > 0) {
                int exon_len = iso->dons[i] - iso->accs[i-1];
                printf("    Exon between intron %d and %d: Length=%d", 
                       i, i+1, exon_len);
                if (exon_len < 50) {  // Assuming min exon is ~50bp
                    printf(" [WARNING: Very short exon!]");
                }
                printf("\n");
            }
        }
        
        // Check terminal exons
        int first_exon = iso->dons[0] - iso->beg;
        int last_exon = iso->end - iso->accs[iso->n_introns-1];
        printf("  First exon length: %d\n", first_exon);
        printf("  Last exon length: %d\n", last_exon);
    }
    printf("=======================\n");
}

/* --------------- Viterbi Algorithm --------------- */

void single_viterbi_algo(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                        Vitbi_algo *vit, Lambda *l, Locus *loc) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    
    // Reset and initialize Viterbi arrays
    reset_viterbi(vit, info);
    initialize_viterbi_from_posterior(vit, pos, info);
    
    int *path           = calloc(info->T, sizeof(int));
    double *path_val    = calloc(info->T, sizeof(double));
    int first_dons      = pos->dons_bps[0];

    for (int t = FLANK; t <= first_dons; t++) {
        path[t] = 0;
        path_val[t] = vit->v[0][t];
    }
    
    double exon_val, intron_val;
    
    for (int t = first_dons + 1; t < info->T - FLANK; t++) {
        exon_val    = vit->v[0][t];
        intron_val  = vit->v[1][t];
        
        if (pos->xi[t][0] != 0.0) {
            vit->v[0][t] = log_sum_sub(exon_val, 0.0, pos->xi[t][0]);
            vit->v[1][t] = log_sum_sub(intron_val, pos->xi[t][0], 0.0);
        }
        else if (pos->xi[t][1] != 0.0) {
            vit->v[0][t] = log_sum_sub(exon_val, pos->xi[t][1], 0.0);
            vit->v[1][t] = log_sum_sub(intron_val, 0.0, pos->xi[t][1]);
        }
        
        if (vit->v[0][t] > vit->v[1][t]) {
            path[t] = 0;
            path_val[t] = vit->v[0][t];
        } else {
            path[t] = 1;
            path_val[t] = vit->v[1][t];
        }
    }
    
    double total_val = 0.0;
    for (int t = FLANK; t < info->T - FLANK; t++) {
        total_val += path_val[t];
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, info->T - FLANK - 1);
        extract_isoform_from_path(path, info, iso);
        iso->val = total_val;
        
        // Validate before adding
        if (validate_isoform(iso, ed)) {
            loc->isoforms[loc->n_isoforms++] = iso;
        } else {
            if (DEBUG) {
                printf("Invalid isoform rejected in single_viterbi_algo\n");
                debug_isoform(iso, "Rejected");
            }
            free_isoform(iso);
        }
    }
    
    free(path);
    free(path_val);
}

void path_restricted_viterbi(Pos_prob *pos, Observed_events *info, Explicit_duration *ed, 
                             Vitbi_algo *vit, Lambda *l, Locus *loc) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    
    // Reset and initialize Viterbi arrays
    reset_viterbi(vit, info);
    initialize_viterbi_from_posterior(vit, pos, info);
    
    int *path               = calloc(info->T, sizeof(int));
    int *last_transition    = calloc(info->T, sizeof(int));
    double *path_val        = calloc(info->T, sizeof(double));
    int first_dons          = pos->dons_bps[0];
    
    for (int t = FLANK; t <= first_dons; t++) {
        path[t] = 0;
        last_transition[t] = FLANK;
        path_val[t] = vit->v[0][t];
    }
    
    double exon_val, intron_val;
    int end_boundary = info->T - FLANK;
    
    for (int t = first_dons + 1; t < end_boundary; t++) {
        exon_val        = vit->v[0][t];
        intron_val      = vit->v[1][t];
        
        int prev_state      = path[t-1];
        int state_duration  = t - last_transition[t-1];
        int remaining_bps   = end_boundary - t;
        
        if (pos->xi[t][0] != 0.0 && prev_state == 0) {
            if (state_duration >= ed->min_len_exon && 
                remaining_bps >= (ed->min_len_intron + ed->min_len_exon)) {
                
                vit->v[0][t] = log_sum_sub(exon_val, 0.0, pos->xi[t][0]);
                vit->v[1][t] = log_sum_sub(intron_val, pos->xi[t][0], 0.0);
            }
        }
        else if (pos->xi[t][1] != 0.0 && prev_state == 1) {
            if (state_duration >= ed->min_len_intron && 
                remaining_bps >= ed->min_len_exon) {
                
                vit->v[0][t] = log_sum_sub(exon_val, pos->xi[t][1], 0.0);
                vit->v[1][t] = log_sum_sub(intron_val, 0.0, pos->xi[t][1]);
            }
        }
        
        if (vit->v[0][t] > vit->v[1][t]) {
            path[t] = 0;
            path_val[t] = vit->v[0][t];
            if (prev_state == 1) {
                last_transition[t] = t;
            } else {
                last_transition[t] = last_transition[t-1];
            }
        } else {
            if (remaining_bps < ed->min_len_exon && prev_state == 1) {
                path[t] = 0;
                path_val[t] = vit->v[0][t];
                last_transition[t] = t;
            } else {
                path[t] = 1;
                path_val[t] = vit->v[1][t];
                if (prev_state == 0) {
                    last_transition[t] = t;
                } else {
                    last_transition[t] = last_transition[t-1];
                }
            }
        }
    }
    
    if (path[end_boundary - 1] == 1) {
        int last_donor_pos = -1;
        for (int t = end_boundary - 1; t >= FLANK; t--) {
            if (path[t] == 1 && (t == FLANK || path[t-1] == 0)) {
                last_donor_pos = t;
                break;
            }
        }
        
        if (last_donor_pos > 0) {
            int forced_acceptor = last_donor_pos + ed->min_len_intron;
            if (forced_acceptor < end_boundary) {
                for (int t = forced_acceptor; t < end_boundary; t++) {
                    path[t] = 0;
                }
            } else {
                for (int t = last_donor_pos; t < end_boundary; t++) {
                    path[t] = 0;
                }
            }
        }
    }
    
    double total_val = 0.0;
    for (int t = FLANK; t < end_boundary; t++) {
        total_val += path_val[t];
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, end_boundary - 1);
        extract_isoform_from_path(path, info, iso);
        iso->val = total_val;
        
        // Validate before adding
        if (validate_isoform(iso, ed)) {
            loc->isoforms[loc->n_isoforms++] = iso;
        } else {
            if (DEBUG) {
                printf("Invalid isoform rejected in path_restricted_viterbi\n");
                debug_isoform(iso, "Rejected");
            }
            free_isoform(iso);
        }
    }
    
    free(path);
    free(last_transition);
    free(path_val);
}

void extract_isoform_from_path(int *path, Observed_events *info, Isoform *iso) {
    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;

    int *temp_dons = malloc(info->T * sizeof(int));
    int *temp_accs = malloc(info->T * sizeof(int));
    int n_dons = 0;
    int n_accs = 0;
    
    for (int t = FLANK + 1; t < info->T - FLANK - 1; t++) {
        if (path[t-1] == 0 && path[t] == 1) {
            temp_dons[n_dons++] = t;
        }
        else if (path[t-1] == 1 && path[t] == 0) {
            temp_accs[n_accs++] = t;
        }
    }
    
    // Set the number of introns
    iso->n_introns = n_dons;

    // BUG FIX: Only allocate and fill arrays if there are introns
    if (n_dons > 0) {
        iso->dons = malloc(n_dons * sizeof(int));
        iso->accs = malloc(n_accs * sizeof(int));  // Note: n_accs should equal n_dons
        
        for (int i = 0; i < n_dons; i++) {
            iso->dons[i] = temp_dons[i];
        }
        for (int i = 0; i < n_accs; i++) {
            iso->accs[i] = temp_accs[i];
        }
    } else {
        // No introns found - set pointers to NULL
        iso->dons = NULL;
        iso->accs = NULL;
    }
    
    // Add debug output
    if (DEBUG && iso->n_introns > 0) {
        debug_isoform(iso, "Just extracted");
    }
    
    free(temp_dons);
    free(temp_accs);
}

/* --------------- Memory Management --------------- */

void free_splice_sites(Pos_prob *pos) {
    free(pos->dons_val);
    free(pos->dons_bps);
    free(pos->accs_val);
    free(pos->accs_bps);
}