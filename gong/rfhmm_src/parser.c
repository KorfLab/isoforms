#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

/* --------------- Auxiliary Function --------------- */

static int count_lines(const char *filename, int skip_header) {
    FILE *file = fopen(filename, "r");
    if (!file) return 0;
    
    int count = 0;
    char line[256];
    
    while (fgets(line, sizeof(line), file)) {
        if (skip_header && line[0] == '%') continue;
        if (line[0] != '\n' && line[0] != '\r') count++;
    }
    
    fclose(file);
    return count;
}

static int detect_kmer_length(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) return 0;
    
    char line[256];
    char seq[20];
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%' || line[0] == '\n' || line[0] == '\r') continue;
        
        // Try to parse sequence
        if (sscanf(line, "%s", seq) == 1) {
            int len = 0;
            for (int i = 0; seq[i] != '\0'; i++) {
                if (seq[i] == 'A' || seq[i] == 'C' || seq[i] == 'G' || seq[i] == 'T') {
                    len++;
                } else {
                    break;
                }
            }
            fclose(file);
            return len;
        }
    }
    
    fclose(file);
    return 0;
}

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

/* --------------- Computation for Transition Prob --------------- */

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

/* --------------- Compute Transition Probability Matrices --------------- */

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

/* --------------- Parser Functions --------------- */

void read_sequence_file(const char *filename, Observed_events *info) {
    if (DEBUG) printf("Reading sequence from: %s\n", filename);
    
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Cannot open sequence file %s\n", filename);
        exit(1);
    }

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    char *buffer = (char*)malloc(file_size + 1);
    fread(buffer, 1, file_size, file);
    buffer[file_size] = '\0';

    // Allocate sequence buffer
    char *sequence = (char*)malloc(file_size + 1);
    size_t seq_index = 0;

    // Parse each line
    char *line = strtok(buffer, "\n");
    while (line != NULL) {
        if (line[0] == '>') {
            if (DEBUG) printf("  Skipping header: %s\n", line);
            line = strtok(NULL, "\n");
            continue;
        }
        
        // Extract valid DNA characters
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T') {
                sequence[seq_index++] = line[i];
            }
        }
        line = strtok(NULL, "\n");
    }
    
    sequence[seq_index] = '\0';
    sequence = (char*)realloc(sequence, seq_index + 1);
    
    info->original_sequence = sequence;
    info->T = seq_index;
    
    free(buffer);
    fclose(file);
    
    if (DEBUG) printf("  Sequence length: %zu bp\n", seq_index);
}

void donor_parser(Lambda *l, char *filename) {
    if (DEBUG) printf("  Parsing donor PWM: %s...", filename);
    
    l->B.don_kmer_len = count_lines(filename, 1);
    if (l->B.don_kmer_len == 0) {
        fprintf(stderr, "ERROR: Cannot determine donor k-mer length\n");
        exit(1);
    }
    
    // Allocate emission matrix for donors
    l->B.dons = malloc(l->B.don_kmer_len * sizeof(double*));
    for (int i = 0; i < l->B.don_kmer_len; i++) {
        l->B.dons[i] = calloc(4, sizeof(double));
    }
    
    // Parse PWM file
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "ERROR: Cannot open donor file %s\n", filename);
        exit(1);
    }
    
    char line[256];
    int row = 0;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;
        
        char *token = strtok(line, " \t\n");
        int col = 0;
        
        while (token && col < 4) {
            l->B.dons[row][col] = atof(token);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }
    
    fclose(file);
    if (DEBUG) printf(" k-mer=%d\n", l->B.don_kmer_len);
}

void acceptor_parser(Lambda *l, char *filename) {
    if (DEBUG) printf("  Parsing acceptor PWM: %s...", filename);
    
    l->B.acc_kmer_len = count_lines(filename, 1);
    if (l->B.acc_kmer_len == 0) {
        fprintf(stderr, "ERROR: Cannot determine acceptor k-mer length\n");
        exit(1);
    }
    
    // Allocate emission matrix for acceptors
    l->B.accs = malloc(l->B.acc_kmer_len * sizeof(double*));
    for (int i = 0; i < l->B.acc_kmer_len; i++) {
        l->B.accs[i] = calloc(4, sizeof(double));
    }
    
    // Parse PWM file
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "ERROR: Cannot open acceptor file %s\n", filename);
        exit(1);
    }
    
    char line[256];
    int row = 0;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;
        
        char *token = strtok(line, " \t\n");
        int col = 0;
        
        while (token && col < 4) {
            l->B.accs[row][col] = atof(token);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }
    
    fclose(file);
    if (DEBUG) printf(" k-mer=%d\n", l->B.acc_kmer_len);
}

void exon_intron_parser(Lambda *l, char *filename, int digit) {
    const char *type = (digit == 0) ? "exon" : "intron";
    if (DEBUG) printf("  Parsing %s MM: %s...", type, filename);
    
    int kmer_len = detect_kmer_length(filename);
    if (kmer_len == 0) {
        fprintf(stderr, "ERROR: Cannot detect k-mer length for %s\n", type);
        exit(1);
    }
    
    if (digit == 0) {
        l->B.exon_kmer_len = kmer_len;
        int size = power(4, kmer_len);
        l->B.exon = calloc(size, sizeof(double));
    } else {
        l->B.intron_kmer_len = kmer_len;
        int size = power(4, kmer_len);
        l->B.intron = calloc(size, sizeof(double));
    }
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "ERROR: Cannot open %s file %s\n", type, filename);
        exit(1);
    }
    
    char line[256];
    char seq[20];
    double prob;
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%' || line[0] == '\n') continue;
        
        if (sscanf(line, "%s %lf", seq, &prob) == 2) {
            int index = 0;
            for (int i = 0; i < kmer_len; i++) {
                int base = 0;
                if (seq[i] == 'C') base = 1;
                else if (seq[i] == 'G') base = 2;
                else if (seq[i] == 'T') base = 3;
                index = index * 4 + base;
            }
            
            if (digit == 0) l->B.exon[index] = prob;
            else l->B.intron[index] = prob;
        }
    }
    
    fclose(file);
    if (DEBUG) printf(" k-mer=%d\n", kmer_len);
}

void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit) {
    const char *type = (digit == 0) ? "exon" : "intron";
    if (DEBUG) printf("  Parsing %s length: %s...", type, filename);
    
    int n_lines = count_lines(filename, 1);
    
    if (digit == 0) {
        ed->exon_len = n_lines;
        ed->exon = calloc(n_lines, sizeof(double));
    } else {
        ed->intron_len = n_lines;
        ed->intron = calloc(n_lines, sizeof(double));
    }
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "ERROR: Cannot open %s duration file %s\n", type, filename);
        exit(1);
    }
    
    char line[256];
    int duration = 0;
    
    if (digit == 0) {
        ed->min_len_exon = -1;
        ed->max_len_exon = 0;
    } else {
        ed->min_len_intron = -1;
        ed->max_len_intron = 0;
    }
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;
        
        double prob = atof(line);
        
        if (digit == 0) {
            ed->exon[duration] = prob;
            if (prob > 0.0 && ed->min_len_exon == -1) {
                ed->min_len_exon = duration;
            }
            if (prob > 0.0) {
                ed->max_len_exon = duration;
            }
        } else {
            ed->intron[duration] = prob;
            if (prob > 0.0 && ed->min_len_intron == -1) {
                ed->min_len_intron = duration;
            }
            if (prob > 0.0) {
                ed->max_len_intron = duration;
            }
        }
        duration++;
    }
    
    fclose(file);
    if (DEBUG) printf(" min=%d, max=%d\n", 
                      digit == 0 ? ed->min_len_exon : ed->min_len_intron,
                      digit == 0 ? ed->max_len_exon : ed->max_len_intron);
}

/* --------------- Memory Cleanup --------------- */

void free_lambda(Lambda *l) {
    // Free emission matrix
    if (l->B.dons) {
        for (int i = 0; i < l->B.don_kmer_len; i++) {
            free(l->B.dons[i]);
        }
        free(l->B.dons);
    }
    
    if (l->B.accs) {
        for (int i = 0; i < l->B.acc_kmer_len; i++) {
            free(l->B.accs[i]);
        }
        free(l->B.accs);
    }
    
    if (l->B.exon) free(l->B.exon);
    if (l->B.intron) free(l->B.intron);
    
    // Free transition matrix
    if (l->A.dons) free(l->A.dons);
    if (l->A.accs) free(l->A.accs);
    
    // Free the rest
    if (l->log_values) free(l->log_values);
}

void free_explicit_duration(Explicit_duration *ed) {
    if (ed->exon) free(ed->exon);
    if (ed->intron) free(ed->intron);
}