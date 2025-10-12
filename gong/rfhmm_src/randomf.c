#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "randomf.h"
#include "model.h"

/* --------------- Hash Table Size Function --------------- */
// using 6k±1 property of finding prime number next to 0.8 load size
static bool is_prime(int n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    
    for (int i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) {
            return false;
        }
    }
    return true;
}

static int next_prime(int n) {
    if (n <= 2) return 2;
    if (n <= 3) return 3;    
    if (n % 2 == 0) n++;

    while (!is_prime(n)) {
        n += 2;
    }
    return n;
}

int next_prime_optimized(int n) {
    if (n <= 2) return 2;
    if (n <= 3) return 3;
    if (n <= 5) return 5;
    
    int remainder = n % 6;
    if (remainder <= 1) {
        int candidate = n - remainder + (remainder == 0 ? 5 : -1);
        if (candidate >= n && is_prime(candidate)) return candidate;
        candidate += 2;
        if (is_prime(candidate)) return candidate;
    } else if (remainder <= 5) {
        int candidate = n - remainder + (remainder <= 1 ? 1 : (remainder <= 5 ? 1 : 5));
        if (remainder > 1) candidate = n - remainder + 1;
        if (candidate < n) candidate += 6;
        
        while (!is_prime(candidate)) {
            candidate += (candidate % 6 == 1) ? 4 : 2;
        }
        return candidate;
    }
    return next_prime(n);
}

int compute_hash_table_size(int locus_size) {
    int min_size = (int)ceil((double)locus_size / 0.8);    
    return next_prime_optimized(min_size);
}

/* --------------- Hash Table Functions --------------- */

IsoformHashTable* create_hash_table(int size) {
    IsoformHashTable *table = malloc(sizeof(IsoformHashTable));
    table->size     = size;
    table->count    = 0;
    table->buckets  = calloc(size, sizeof(HashNode*));
    return table;
}

void free_hash_table(IsoformHashTable *table) {
    if (!table) return;
    
    for (int i = 0; i < table->size; i++) {
        HashNode *current = table->buckets[i];
        while (current) {
            HashNode *next = current->next;
            free(current);
            current = next;
        }
    }
    free(table->buckets);
    free(table);
}

unsigned long compute_isoform_hash(Isoform *iso, int table_size) {
    if (iso->n_introns == 0) {
        return 0;
    }
    
    unsigned long hash = 0;
    unsigned long sum_donors = 0;
    unsigned long sum_acceptors = 0;
    
    for (int i = 0; i < iso->n_introns; i++) {
        sum_donors += iso->dons[i];
        sum_acceptors += iso->accs[i];
    }
    
    hash = (iso->n_introns * 31UL +
            iso->n_introns * 47UL +
            sum_donors * 73UL +
            sum_acceptors * 101UL)
            % table_size;
    
    return hash;
}

int isoforms_are_identical(Isoform *iso1, Isoform *iso2) {
    if (iso1->n_introns != iso2->n_introns) {
        return 0;
    }

    if (iso1->n_introns == 0) {
        return 1;
    }

    for (int i = 0; i < iso1->n_introns; i++) {
        if (iso1->dons[i] != iso2->dons[i] || 
            iso1->accs[i] != iso2->accs[i]) {
            return 0;
        }
    }
    return 1;
}

int isoform_exists_in_hash(IsoformHashTable *table, Isoform *new_iso) {
    if (!table || !new_iso) return 0;
    
    unsigned long hash = compute_isoform_hash(new_iso, table->size);
    HashNode *current = table->buckets[hash];    
    while (current != NULL) {
        if (isoforms_are_identical(current->isoform, new_iso)) {
            return 1;
        }
        current = current->next;
    }
    
    return 0;
}

void insert_isoform_to_hash(IsoformHashTable *table, Isoform *iso) {
    if (!table || !iso) return;
    
    unsigned long hash  = compute_isoform_hash(iso, table->size);
    HashNode *new_node  = malloc(sizeof(HashNode));
    new_node->isoform   = iso;
    new_node->next      = table->buckets[hash];
    // add to head
    table->buckets[hash] = new_node;
    table->count++;
}

/* --------------- Initialize Random Forest--------------- */

RandomForest* create_random_forest(Pos_prob *pos, Locus *loc, int node_size, float mtry) {
    RandomForest *rf = malloc(sizeof(RandomForest));
    
    rf->mtry                = mtry;
    rf->sample_size         = pos->dons+pos->accs;
    rf->all_sites           = malloc(rf->sample_size * sizeof(SpliceSite));
    rf->node_size           = node_size;
    rf->hash_table          = create_hash_table(compute_hash_table_size(loc->capacity));
    // Append original dataset
    int idx = 0;
    for (int i = 0; i < pos->dons; i++) {
        rf->all_sites[idx].pos = pos->dons_bps[i];
        rf->all_sites[idx].typ = 0;
        rf->all_sites[idx].val = pos->dons_val[i];
        idx++;
    }
    for (int i = 0; i < pos->accs; i++) {
        rf->all_sites[idx].pos = pos->accs_bps[i];
        rf->all_sites[idx].typ = 1;
        rf->all_sites[idx].val = pos->accs_val[i];
        idx++;
    }
    srand(time(NULL));
    return rf;
}

/* --------------- Bootstrap Sampling --------------- */
// bootstrap sampling (replacement = True)
static SpliceSite* bootstrap_sample(SpliceSite *sites, int n_sites) {
    if (n_sites <= 0 || sites == NULL) {
        return NULL;
    }
    SpliceSite *sample = malloc(n_sites * sizeof(SpliceSite));
    for (int i = 0; i < n_sites; i++) {
        sample[i] = sites[rand() % n_sites];
    }
    return sample;
}
// for mtry sampling
static int remove_duplicates(SpliceSite *sites, int n_sites, SpliceSite **unique_sites) {
    if (n_sites == 0 || sites == NULL) {
        *unique_sites = NULL;
        return 0;
    }
    
    int *seen = calloc(n_sites, sizeof(int));
    int unique_count = 0;
    
    for (int i = 0; i < n_sites; i++) {
        if (seen[i]) continue;
        
        seen[i] = 1;
        unique_count++;
        
        for (int j = i + 1; j < n_sites; j++) {
            if (!seen[j] && 
                sites[i].pos == sites[j].pos && 
                sites[i].typ == sites[j].typ) {
                seen[j] = 1;
            }
        }
    }
    
    *unique_sites = malloc(unique_count * sizeof(SpliceSite));
    int idx = 0;
    
    memset(seen, 0, n_sites * sizeof(int));
    for (int i = 0; i < n_sites; i++) {
        if (seen[i]) continue;
        
        (*unique_sites)[idx++] = sites[i];
        seen[i] = 1;
        
        for (int j = i + 1; j < n_sites; j++) {
            if (!seen[j] && 
                sites[i].pos == sites[j].pos && 
                sites[i].typ == sites[j].typ) {
                seen[j] = 1;
            }
        }
    }
    
    free(seen);
    return unique_count;
}

/* --------------- Sample without replacement --------------- */
// Optimized Fisher-Yates Shuffle
static void fisher_yates_shuffle(int *array, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j       = rand() % (i + 1);
        int temp    = array[i];
        array[i]    = array[j];
        array[j]    = temp;
    }
}

/* --------------- Splitting Criteria --------------- */

static double compute_mse(SpliceSite *sites, int n_sites) {
    if (n_sites == 0) return 0.0;
    
    double mean = 0.0;
    for (int i = 0; i < n_sites; i++) {
        mean += sites[i].val;
    }
    mean /= n_sites;
    
    double mse = 0.0;
    for (int i = 0; i < n_sites; i++) {
        double diff = sites[i].val - mean;
        mse += diff * diff;
    }
    return mse / n_sites;
}

static double compute_gini(SpliceSite *sites, int n_sites) {
    if (n_sites == 0) return 0.0;
    
    int donors = 0, acceptors = 0;
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].typ == 0) donors++;
        else acceptors++;
    }
    
    double p_donor      = (double)donors / n_sites;
    double p_acceptor   = (double)acceptors / n_sites;
    
    return 1.0 - (p_donor * p_donor + p_acceptor * p_acceptor);
}

static int compare_sites_by_val(const void *a, const void *b) {
    SpliceSite *sa = (SpliceSite*)a;
    SpliceSite *sb = (SpliceSite*)b;
    if      (sa->val < sb->val) return -1;
    else if (sa->val > sb->val) return 1;
    else return 0;
}

/* --------------- Recursive Splitting --------------- */

static int find_best_split(SpliceSite *sites, int n_sites, double *best_threshold,
                          int node_size, float mtry) {

    if (n_sites < node_size * 2 || sites == NULL) {
        return 0;
    }
    
    int subset_size = (int)(n_sites * mtry);
    if (subset_size < 1) subset_size = 1;
    if (subset_size > n_sites) subset_size = n_sites;
    
    int *indices = malloc(n_sites * sizeof(int));
    for (int i = 0; i < n_sites; i++) {
        indices[i] = i;
    }
    
    fisher_yates_shuffle(indices, n_sites);    
    SpliceSite *subset = malloc(subset_size * sizeof(SpliceSite));
    for (int i = 0; i < subset_size; i++) {
        subset[i] = sites[indices[i]];
    }
    free(indices);
    
    qsort(subset, subset_size, sizeof(SpliceSite), compare_sites_by_val);
    
    double parent_mse   = compute_mse(subset, subset_size);
    double max_gain     = -1.0;
    *best_threshold     = 0.0;
    int found_split     = 0;
    
    for (int i = 1; i < subset_size; i++) {
        // Skip if same value as previous (no split possible)
        if (subset[i].val == subset[i-1].val) continue;
        
        double threshold = (subset[i-1].val + subset[i].val) / 2.0;
        
        // Count how split would divide the FULL node dataset
        int left_count = 0, right_count = 0;
        for (int j = 0; j < n_sites; j++) {
            if (sites[j].val < threshold) left_count++;
            else right_count++;
        }
        
        // Enforce minimum node size
        if (left_count < node_size || right_count < node_size) continue;

        // Calculate MSE gain
        double left_mse     = compute_mse(subset, i);
        double right_mse    = compute_mse(&subset[i], subset_size - i);
        double gain         = parent_mse - left_mse - right_mse;

        SpliceSite *temp_left   = malloc(left_count  * sizeof(SpliceSite));
        SpliceSite *temp_right  = malloc(right_count * sizeof(SpliceSite));
        
        int l_idx = 0, r_idx = 0;
        for (int j = 0; j < n_sites; j++) {
            if (sites[j].val < threshold) {
                temp_left[l_idx++]  = sites[j];
            } else {
                temp_right[r_idx++] = sites[j];
            }
        }
        
        // Using Gini impurity to avoid pure dons and accs
        double left_gini    = compute_gini(temp_left,  left_count);
        double right_gini   = compute_gini(temp_right, right_count);
        
        free(temp_left);
        free(temp_right);
        
        // Skip if either child would be pure (all donors or all acceptors)
        if (left_gini == 0.0 || right_gini == 0.0) continue;

        if (gain > max_gain) {
            max_gain = gain;
            *best_threshold = threshold;
            found_split = 1;
        }
    }
    
    free(subset);
    return found_split;
}

/* --------------- Viterbi Algorithm --------------- */

void viterbi_on_subset(SpliceSite *sites, int n_sites, Observed_events *info,
                      Explicit_duration *ed, Lambda *l, Locus *loc, 
                      Vitbi_algo *vit, int node_size,
                      IsoformHashTable *hash_table) {
    
    if (n_sites <= node_size || sites == NULL) {
        return;
    }
    
    Pos_prob subset_pos;
    int n_dons = 0, n_accs = 0;
    
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].typ == 0) n_dons++;
        else n_accs++;
    }
    
    subset_pos.dons_bps = malloc(n_dons * sizeof(int));
    subset_pos.dons_val = malloc(n_dons * sizeof(double));
    subset_pos.accs_bps = malloc(n_accs * sizeof(int));
    subset_pos.accs_val = malloc(n_accs * sizeof(double));
    
    int d_idx = 0, a_idx = 0;
    for (int i = 0; i < n_sites; i++) {
        if (sites[i].typ == 0) {
            subset_pos.dons_bps[d_idx] = sites[i].pos;
            subset_pos.dons_val[d_idx] = sites[i].val;
            d_idx++;
        } else {
            subset_pos.accs_bps[a_idx] = sites[i].pos;
            subset_pos.accs_val[a_idx] = sites[i].val;
            a_idx++;
        }
    }
    
    subset_pos.dons = n_dons;
    subset_pos.accs = n_accs;
    
    int T = info->T;
    subset_pos.xi = malloc(T * sizeof(double*));
    for (int i = 0; i < T; i++) {
        subset_pos.xi[i] = calloc(2, sizeof(double));
    }
    
    for (int i = 0; i < n_dons; i++) {
        subset_pos.xi[subset_pos.dons_bps[i]][0] = subset_pos.dons_val[i];
    }
    for (int i = 0; i < n_accs; i++) {
        subset_pos.xi[subset_pos.accs_bps[i]][1] = subset_pos.accs_val[i];
    }
    // create temp viterbi path
    reset_viterbi(vit, info);    
    initialize_viterbi_from_posterior(vit, &subset_pos, info);
    
    int prev_count = loc->n_isoforms;
    
    path_restricted_viterbi(&subset_pos, info, ed, vit, l, loc);
    // check duplicate isoform
    if (loc->n_isoforms > prev_count) {
        Isoform *new_iso = loc->isoforms[loc->n_isoforms - 1];
        // insert into hash table
        if (isoform_exists_in_hash(hash_table, new_iso)) {
            free_isoform(new_iso);
            loc->n_isoforms--;
        } else {
            insert_isoform_to_hash(hash_table, new_iso);
        }
    }
    // clean up
    free(subset_pos.dons_bps);
    free(subset_pos.dons_val);
    free(subset_pos.accs_bps);
    free(subset_pos.accs_val);
    for (int i = 0; i < T; i++) {
        free(subset_pos.xi[i]);
    }
    free(subset_pos.xi);
}

/* --------------- Build Tree with Direct Viterbi --------------- */

void build_tree_with_viterbi(SpliceSite *sites, int n_sites, RandomForest *rf,
                             Observed_events *info, Explicit_duration *ed, 
                             Lambda *l, Locus *loc, Vitbi_algo *vit) {
    
    SpliceSite *unique_sites = NULL;
    int unique_count = remove_duplicates(sites, n_sites, &unique_sites);
    
    if (DEBUG && unique_count < n_sites) {
        printf("  Node had %d sites, reduced to %d unique sites\n", 
               n_sites, unique_count);
    }
    
    if (unique_count >= rf->node_size) {
        viterbi_on_subset(unique_sites, unique_count, info, ed, l, loc, vit,
                        rf->node_size, rf->hash_table);
    }
    
    double threshold;
    if (!find_best_split(unique_sites, unique_count, &threshold, rf->node_size, rf->mtry)) {
        free(unique_sites);
        return;
    }
    
    // Split data based on threshold
    SpliceSite *left_sites  = malloc(unique_count * sizeof(SpliceSite));
    SpliceSite *right_sites = malloc(unique_count * sizeof(SpliceSite));
    int left_count = 0, right_count = 0;
    
    for (int i = 0; i < unique_count; i++) {
        if (unique_sites[i].val < threshold) {
            left_sites[left_count++] = unique_sites[i];
        } else {
            right_sites[right_count++] = unique_sites[i];
        }
    }
    
    // Recursion on children
    if (left_count > 0) {
        build_tree_with_viterbi(left_sites, left_count, rf, info, ed, l,
                                loc, vit);
    }
    if (right_count > 0) {
        build_tree_with_viterbi(right_sites, right_count, rf, info, ed, l,
                                loc, vit);
    }
    
    free(unique_sites);
    free(left_sites);
    free(right_sites);
}

/* --------------- Execute Isoform Forest --------------- */

void generate_isoforms_random_forest(RandomForest *rf, Observed_events *info,
                                     Explicit_duration *ed, Lambda *l, 
                                     Locus *loc, Vitbi_algo *vit) {

    int trees_without_new_isoforms  = 0;
    int max_trees_without_progress  = 100;
    int tree_count                  = 0;
    
    while (loc->n_isoforms < loc->capacity) {
        tree_count++;
        int prev_isoform_count = loc->n_isoforms;
        
        // Bootstrap sample
        SpliceSite *bootstrap = bootstrap_sample(rf->all_sites, rf->sample_size);
        
        if (DEBUG && tree_count % 10 == 0) {
            printf("Tree %d: %d unique isoforms found so far\n", 
                   tree_count, loc->n_isoforms);
        }
        
        // Build tree and collect isoforms
        build_tree_with_viterbi(bootstrap, rf->sample_size, rf, info, ed, l, 
                               loc, vit);
        free(bootstrap);
        
        // Check if we found new isoforms
        if (loc->n_isoforms == prev_isoform_count) {
            trees_without_new_isoforms++;
            if (trees_without_new_isoforms >= max_trees_without_progress) {
                if (DEBUG) {
                    printf("No new isoforms found in last %d trees. Stopping early.\n", 
                           max_trees_without_progress);
                    printf("Total trees generated: %d\n", tree_count);
                }
                break;
            }
        } else {
            trees_without_new_isoforms = 0;  // Reset counter
        }
        
        // Check if capacity reached
        if (loc->n_isoforms >= loc->capacity) {
            if (DEBUG) {
                printf("Reached isoform capacity (%d) after %d trees\n", 
                       loc->capacity, tree_count);
            }
            break;
        }
    }
    
    if (DEBUG) {
        printf("Random Forest generation complete:\n");
        printf("  Trees generated: %d\n", tree_count);
        printf("  Unique isoforms found: %d\n", loc->n_isoforms);
    }
}

/* --------------- Cleanup --------------- */

void free_random_forest(RandomForest *rf) {
    if (rf) {
        if (rf->all_sites) {
            free(rf->all_sites);
        }
        if (rf->hash_table) {
            if (DEBUG) {
                print_hash_table_stats(rf->hash_table);
            }
            free_hash_table(rf->hash_table);
        }
        free(rf);
    }
}

void print_hash_table_stats(IsoformHashTable *table) {
    if (!table) return;
    
    int used_buckets = 0;
    int max_chain_length = 0;
    int total_chain_length = 0;
    
    for (int i = 0; i < table->size; i++) {
        if (table->buckets[i] != NULL) {
            used_buckets++;
            int chain_length = 0;
            HashNode *current = table->buckets[i];
            while (current) {
                chain_length++;
                current = current->next;
            }
            total_chain_length += chain_length;
            if (chain_length > max_chain_length) {
                max_chain_length = chain_length;
            }
        }
    }
    
    printf("Hash Table Statistics:\n");
    printf("  Table size: %d\n", table->size);
    printf("  Total items: %d\n", table->count);
    printf("  Load factor: %.2f\n", (double)table->count / table->size);
    printf("  Used buckets: %d (%.1f%%)\n", used_buckets, 
           100.0 * used_buckets / table->size);
    printf("  Max chain length: %d\n", max_chain_length);
    printf("  Avg chain length: %.2f\n", 
           used_buckets > 0 ? (double)total_chain_length / used_buckets : 0.0);
}