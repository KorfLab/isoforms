#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"

/* --------------- Auxiliary Functions --------------- */

double total_prob(double *array, int length) {
    double value = 0.0;
    for (int i = 0 ; i < length ; i ++) {
        if (array[i] == 0.0) {
            value = 0.0;
            break;
        }
        value += log(array[i]);
    }

    if (value != 0.0)   value = exp(value);
    return value;
}

/* --------------- Log-Space Transition Computation --------------- */

static void compute_log_transition_recursive_don(Lambda *l, int depth, double log_sum) {
    if (depth == l->B.don_kmer_len) {
        int idx        = base4_to_int(l->A.pos, 0, l->B.don_kmer_len);
        l->A.dons[idx] = log_sum;
        return;
    }

    for (int i = 0; i < 4; i++) {
        l->A.pos[depth] = i;
        double log_prob = l->B.dons[depth][i];
        // In log-space: sum instead of multiply
        double new_sum  = (log_prob <= -99.0) ? LOG_ZERO : (log_sum + log_prob);
        compute_log_transition_recursive_don(l, depth + 1, new_sum);
    }
}

static void compute_log_transition_recursive_acc(Lambda *l, int depth, double log_sum) {
    if (depth == l->B.acc_kmer_len) {
        int idx        = base4_to_int(l->A.pos, 0, l->B.acc_kmer_len);
        l->A.accs[idx] = log_sum;
        return;
    }

    for (int i = 0; i < 4; i++) {
        l->A.pos[depth] = i;
        double log_prob = l->B.accs[depth][i];
        double new_sum  = (log_prob <= -99.0) ? LOG_ZERO : (log_sum + log_prob);
        compute_log_transition_recursive_acc(l, depth + 1, new_sum);
    }
}

/* --------------- Probability-Space Transition Computation --------------- */

static void initialize_donor_transition_matrix_recursive(Lambda *l, int depth) {
    if (depth == l->B.don_kmer_len) {
        int     idx     = base4_to_int(l->A.pos, 0, l->B.don_kmer_len);
        double  val     = total_prob(l->A.prob, l->B.don_kmer_len);
        l->A.dons[idx]  = val;
        return;
    }

    for (int i = 0; i < 4 ; i++) {
        double prob         = l->B.dons[depth][i];
        l->A.prob[depth]    = prob;
        l->A.pos[depth]     = i;
        initialize_donor_transition_matrix_recursive(l, depth+1);
    }
}

static void initialize_acceptor_transition_matrix_recursive(Lambda *l, int depth) {
    if (depth == l->B.acc_kmer_len) {
        int     idx     = base4_to_int(l->A.pos, 0, l->B.acc_kmer_len);
        double  val     = total_prob(l->A.prob, l->B.acc_kmer_len);
        l->A.accs[idx]  = val;
        return;
    }

    for (int i = 0; i < 4 ; i++) {
        double prob         = l->B.accs[depth][i];
        l->A.prob[depth]    = prob;
        l->A.pos[depth]     = i;
        initialize_acceptor_transition_matrix_recursive(l, depth+1);
    }
}

/* --------------- Public Functions --------------- */

void compute_transition_matrices(Lambda *l) {
    if (DEBUG) printf("Computing transition probability matrices...\n");

    // Allocate transition matrices
    l->A.don_size   = power(4, l->B.don_kmer_len);
    l->A.acc_size   = power(4, l->B.acc_kmer_len);
    l->A.dons       = calloc(l->A.don_size, sizeof(double));
    l->A.accs       = calloc(l->A.acc_size, sizeof(double));

    // Temporary arrays for recursive computation
    l->A.pos        = calloc(l->B.acc_kmer_len, sizeof(int));  // Use larger of the two
    l->A.prob       = calloc(l->B.acc_kmer_len, sizeof(double));

    // Compute donor transition matrix
    if (DEBUG) printf("  Computing donor transition matrix...");
    initialize_donor_transition_matrix_recursive(l, 0);
    if (DEBUG) printf(" Done\n");

    // Compute acceptor transition matrix
    if (DEBUG) printf("  Computing acceptor transition matrix...");
    initialize_acceptor_transition_matrix_recursive(l, 0);
    if (DEBUG) printf(" Done\n");

    // Clean up temporary arrays
    free(l->A.pos);
    free(l->A.prob);
    l->A.pos = NULL;
    l->A.prob = NULL;

    if (DEBUG) {
        // Verify non-zero entries
        int non_zero_dons = 0, non_zero_accs = 0;
        for (int i = 0; i < l->A.don_size; i++) {
            if (l->A.dons[i] > 0.0) non_zero_dons++;
        }
        for (int i = 0; i < l->A.acc_size; i++) {
            if (l->A.accs[i] > 0.0) non_zero_accs++;
        }
        printf("  Donor matrix: %d non-zero entries out of %d\n", non_zero_dons, l->A.don_size);
        printf("  Acceptor matrix: %d non-zero entries out of %d\n", non_zero_accs, l->A.acc_size);
    }
}

void compute_transition_matrices_log(Lambda *l) {
    if (DEBUG) printf("Computing log-space transition matrices...\n");

    // Allocate transition matrices
    l->A.don_size = power(4, l->B.don_kmer_len);
    l->A.acc_size = power(4, l->B.acc_kmer_len);
    l->A.dons     = malloc(l->A.don_size * sizeof(double));
    l->A.accs     = malloc(l->A.acc_size * sizeof(double));

    // Initialize with LOG_ZERO
    for (int i = 0; i < l->A.don_size; i++) l->A.dons[i] = LOG_ZERO;
    for (int i = 0; i < l->A.acc_size; i++) l->A.accs[i] = LOG_ZERO;

    // Temporary position array for recursive computation
    int max_kmer = (l->B.don_kmer_len > l->B.acc_kmer_len) ?
                    l->B.don_kmer_len : l->B.acc_kmer_len;
    l->A.pos  = calloc(max_kmer, sizeof(int));
    l->A.prob = NULL;

    // Compute donor transition matrix in log-space
    if (DEBUG) printf("  Computing donor transition matrix...");
    compute_log_transition_recursive_don(l, 0, 0.0);
    if (DEBUG) printf(" Done\n");

    // Compute acceptor transition matrix in log-space
    if (DEBUG) printf("  Computing acceptor transition matrix...");
    compute_log_transition_recursive_acc(l, 0, 0.0);
    if (DEBUG) printf(" Done\n");

    // Clean up
    free(l->A.pos);
    l->A.pos = NULL;

    if (DEBUG) {
        int non_zero_dons = 0, non_zero_accs = 0;
        for (int i = 0; i < l->A.don_size; i++) {
            if (!IS_LOG_ZERO(l->A.dons[i])) non_zero_dons++;
        }
        for (int i = 0; i < l->A.acc_size; i++) {
            if (!IS_LOG_ZERO(l->A.accs[i])) non_zero_accs++;
        }
        printf("  Donor matrix: %d non-zero entries out of %d\n",
               non_zero_dons, l->A.don_size);
        printf("  Acceptor matrix: %d non-zero entries out of %d\n",
               non_zero_accs, l->A.acc_size);
    }
}
