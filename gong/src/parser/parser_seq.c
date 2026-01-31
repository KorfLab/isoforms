#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"

/* --------------- Sequence File Parsing --------------- */

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
