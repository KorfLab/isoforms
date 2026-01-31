#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

/* --------------- Auxiliary Functions --------------- */

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

/* --------------- Markov Model Parsing --------------- */

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
