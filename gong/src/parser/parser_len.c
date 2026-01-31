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

/* --------------- Length Distribution Parsing --------------- */

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
