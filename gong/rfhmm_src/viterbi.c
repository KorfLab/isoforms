#include <limits.h>
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

/* --------------- Viterbi State Management --------------- */

void allocate_vit(Vitbi_algo *vit, Observed_events *info) {
    if (!vit || !info || info->T <= 0) {
        return;
    }

    vit->v = malloc(HS * sizeof(double*));

    if (!vit->v) {
        fprintf(stderr, "Failed to allocate Viterbi state arrays.\n");
        exit(EXIT_FAILURE);
    }

    for (int state = 0; state < HS; state++) {
        vit->v[state] = calloc(info->T, sizeof(double));
        if (!vit->v[state]) {
            fprintf(stderr, "Failed to allocate Viterbi state rows.\n");
            exit(EXIT_FAILURE);
        }
    }

    vit->exon = 0.0;
    vit->intron = 0.0;
}

void free_vit(Vitbi_algo *vit, Observed_events *info) {
    (void)info;

    if (!vit || !vit->v) {
        return;
    }

    for (int state = 0; state < HS; state++) {
        free(vit->v[state]);
        vit->v[state] = NULL;
    }
    free(vit->v);
    vit->v = NULL;

    vit->exon = 0.0;
    vit->intron = 0.0;
}

void reset_viterbi(Vitbi_algo *vit, Observed_events *info) {
    if (!vit || !vit->v || !info) {
        return;
    }

    for (int state = 0; state < HS; state++) {
        memset(vit->v[state], 0, sizeof(double) * info->T);
    }

    vit->exon = 0.0;
    vit->intron = 0.0;
}

void initialize_viterbi_from_posterior(Vitbi_algo *vit, Pos_prob *pos, Observed_events *info) {
    if (!vit || !vit->v || !pos || !pos->xi || !info) {
        return;
    }

    for (int t = 0; t < info->T; t++) {
        double *posterior_states = pos->xi[t];
        for (int state = 0; state < HS; state++) {
            double value = 0.0;
            if (posterior_states) {
                value = posterior_states[state];
            }
            vit->v[state][t] = value;
        }
    }
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
    
    int end_boundary = info->T - FLANK;

    if (path[end_boundary - 1] == 1) {
        int last_donor_pos = -1;
        for (int t = end_boundary - 1; t >= FLANK; t--) {
            if (path[t] == 1 && (t == FLANK || path[t-1] == 0)) {
                last_donor_pos = t;
                break;
            }
        }

        if (last_donor_pos >= 0) {
            int forced_acceptor = last_donor_pos + ed->min_len_intron;
            int start_zero = forced_acceptor < end_boundary ? forced_acceptor : end_boundary;
            for (int t = start_zero; t < end_boundary; t++) {
                path[t] = 0;
            }
        }
    }

    double total_val = 0.0;
    for (int t = FLANK; t < end_boundary; t++) {
        total_val += path_val[t];
    }
    
    if (loc->n_isoforms < loc->capacity) {
        Isoform *iso = create_isoform(FLANK, end_boundary-1);
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

typedef struct {
    int pos;
    double val;
} WeightedSite;

static int compare_site_positions(const void *a, const void *b) {
    const WeightedSite *sa = (const WeightedSite*)a;
    const WeightedSite *sb = (const WeightedSite*)b;
    if (sa->pos < sb->pos) return -1;
    if (sa->pos > sb->pos) return 1;
    return 0;
}

static void assemble_restricted_isoform(Pos_prob *pos, Observed_events *info,
                                        Explicit_duration *ed, Isoform *iso) {
    (void)info;

    int min_exon   = (ed->min_len_exon   < 0) ? 0 : ed->min_len_exon;
    int min_intron = (ed->min_len_intron < 0) ? 0 : ed->min_len_intron;
    int max_exon   = (ed->max_len_exon   > 0 && ed->max_len_exon >= min_exon)
                     ? ed->max_len_exon : INT_MAX;
    int max_intron = (ed->max_len_intron > 0 && ed->max_len_intron >= min_intron)
                     ? ed->max_len_intron : INT_MAX;

    int donor_count = pos->dons;
    int acceptor_count = pos->accs;
    if (donor_count <= 0 || acceptor_count <= 0) {
        iso->n_introns = 0;
        iso->val = 0.0;
        return;
    }

    WeightedSite *donors = malloc(donor_count * sizeof(WeightedSite));
    WeightedSite *acceptors = malloc(acceptor_count * sizeof(WeightedSite));
    if (!donors || !acceptors) {
        free(donors);
        free(acceptors);
        iso->n_introns = 0;
        iso->val = 0.0;
        return;
    }

    for (int i = 0; i < donor_count; i++) {
        donors[i].pos = pos->dons_bps[i];
        donors[i].val = pos->dons_val ? pos->dons_val[i] : 0.0;
    }
    for (int i = 0; i < acceptor_count; i++) {
        acceptors[i].pos = pos->accs_bps[i];
        acceptors[i].val = pos->accs_val ? pos->accs_val[i] : 0.0;
    }

    qsort(donors, donor_count, sizeof(WeightedSite), compare_site_positions);
    qsort(acceptors, acceptor_count, sizeof(WeightedSite), compare_site_positions);

    int capacity = (donor_count < acceptor_count) ? donor_count : acceptor_count;
    int *selected_dons = malloc(capacity * sizeof(int));
    int *selected_accs = malloc(capacity * sizeof(int));
    double *selected_don_vals = malloc(capacity * sizeof(double));
    double *selected_acc_vals = malloc(capacity * sizeof(double));

    if (!selected_dons || !selected_accs || !selected_don_vals || !selected_acc_vals) {
        free(selected_dons);
        free(selected_accs);
        free(selected_don_vals);
        free(selected_acc_vals);
        free(donors);
        free(acceptors);
        iso->n_introns = 0;
        iso->val = 0.0;
        return;
    }

    int pair_count = 0;
    double total_score = 0.0;
    int prev_acceptor = iso->beg;
    int acc_index = 0;

    for (int d = 0; d < donor_count; d++) {
        int donor_pos = donors[d].pos;
        if (donor_pos >= iso->end) {
            break;
        }

        if (donor_pos <= prev_acceptor) {
            continue;
        }

        int exon_len = donor_pos - prev_acceptor;
        if (exon_len < min_exon) {
            continue;
        }
        if (exon_len > max_exon) {
            break;
        }

        while (acc_index < acceptor_count && acceptors[acc_index].pos <= donor_pos) {
            acc_index++;
        }

        int chosen_idx = -1;
        for (int a = acc_index; a < acceptor_count; a++) {
            int acc_pos = acceptors[a].pos;
            if (acc_pos > iso->end) {
                break;
            }
            int intron_len = acc_pos - donor_pos;
            if (intron_len < min_intron) {
                continue;
            }
    
            if (intron_len > max_intron) {
                break;
            }

            chosen_idx = a;
            break;
        }

        if (chosen_idx == -1) {
            continue;
        }

        selected_dons[pair_count] = donor_pos;
        selected_accs[pair_count] = acceptors[chosen_idx].pos;
        selected_don_vals[pair_count] = donors[d].val;
        selected_acc_vals[pair_count] = acceptors[chosen_idx].val;
        total_score += donors[d].val + acceptors[chosen_idx].val;

        prev_acceptor = acceptors[chosen_idx].pos;
        acc_index = chosen_idx + 1;
        pair_count++;

        if (pair_count == capacity) {
            break;
        }
    }

    if (pair_count > 0) {
        int final_exon = iso->end - selected_accs[pair_count - 1];
        while (pair_count > 0 && final_exon < min_exon) {
            pair_count--;
            total_score -= selected_don_vals[pair_count] + selected_acc_vals[pair_count];
            final_exon = (pair_count > 0)
                         ? iso->end - selected_accs[pair_count - 1]
                         : iso->end - iso->beg;
        }
    }

    if (pair_count == 0) {
        iso->n_introns = 0;
        iso->val = 0.0;
        free(selected_dons);
        free(selected_accs);
        free(selected_don_vals);
        free(selected_acc_vals);
        free(donors);
        free(acceptors);
        return;
    }

    iso->dons = malloc(pair_count * sizeof(int));
    iso->accs = malloc(pair_count * sizeof(int));
    if (!iso->dons || !iso->accs) {
        free(iso->dons);
        free(iso->accs);
        iso->dons = NULL;
        iso->accs = NULL;
        iso->n_introns = 0;
        iso->val = 0.0;
        free(selected_dons);
        free(selected_accs);
        free(selected_don_vals);
        free(selected_acc_vals);
        free(donors);
        free(acceptors);
        return;
    }

    memcpy(iso->dons, selected_dons, pair_count * sizeof(int));
    memcpy(iso->accs, selected_accs, pair_count * sizeof(int));
    iso->n_introns = pair_count;
    iso->val = total_score;

    free(selected_dons);
    free(selected_accs);
    free(selected_don_vals);
    free(selected_acc_vals);
    free(donors);
    free(acceptors);
}

void path_restricted_viterbi(Pos_prob *pos, Observed_events *info, Explicit_duration *ed,
                             Vitbi_algo *vit, Lambda *l, Locus *loc) {
    (void)vit;
    (void)l;

    if (!pos || !info || !ed || !loc) {
        return;
    }

    if (loc->n_isoforms >= loc->capacity) {
        return;
    }

    int FLANK = (info->flank != 0) ? info->flank : DEFAULT_FLANK;
    int end_boundary = info->T - FLANK;

    Isoform *iso = create_isoform(FLANK, end_boundary - 1);
    assemble_restricted_isoform(pos, info, ed, iso);

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
    
    int max_pairs = (n_dons < n_accs) ? n_dons : n_accs;
    if (max_pairs == 0) {
        if (DEBUG) {
            printf("Extracted transitions from path: donors=%d, acceptors=%d, paired=0\n",
                   n_dons, n_accs);
        }
        iso->n_introns = 0;
        iso->dons = NULL;
        iso->accs = NULL;
        free(temp_dons);
        free(temp_accs);
        return;
    }

    iso->dons = malloc(max_pairs * sizeof(int));
    iso->accs = malloc(max_pairs * sizeof(int));
    if (!iso->dons || !iso->accs) {
        free(iso->dons);
        free(iso->accs);
        iso->dons = NULL;
        iso->accs = NULL;
        iso->n_introns = 0;
        free(temp_dons);
        free(temp_accs);
        return;
    }

    int don_idx = 0;
    int acc_idx = 0;
    int pair_idx = 0;

    while (don_idx < n_dons && acc_idx < n_accs) {
        if (temp_dons[don_idx] < temp_accs[acc_idx]) {
            iso->dons[pair_idx] = temp_dons[don_idx];
            iso->accs[pair_idx] = temp_accs[acc_idx];
            pair_idx++;
            don_idx++;
            acc_idx++;
        } else {
            acc_idx++;
        }
    }

    iso->n_introns = pair_idx;

    if (DEBUG) {
        printf("Extracted transitions from path: donors=%d, acceptors=%d, paired=%d\n",
               n_dons, n_accs, iso->n_introns);
    }

    if (iso->n_introns == 0) {
        free(iso->dons);
        free(iso->accs);
        iso->dons = NULL;
        iso->accs = NULL;
    } else if (iso->n_introns < max_pairs) {
        int new_size = iso->n_introns * sizeof(int);
        int *resized_dons = realloc(iso->dons, new_size);
        int *resized_accs = realloc(iso->accs, new_size);
        if (resized_dons) iso->dons = resized_dons;
        if (resized_accs) iso->accs = resized_accs;
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