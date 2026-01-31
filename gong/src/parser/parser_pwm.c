#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

/* --------------- Auxiliary Functions --------------- */

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

/* --------------- PWM Parsing --------------- */

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
