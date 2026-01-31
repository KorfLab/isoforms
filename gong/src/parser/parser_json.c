#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "model.h"
#include "cJSON.h"

/* --------------- Auxiliary Functions --------------- */

static int base_to_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

static int kmer_to_index(const char *kmer, int len) {
    int index = 0;
    for (int i = 0; i < len; i++) {
        int base = base_to_index(kmer[i]);
        if (base < 0) return -1;
        index = index * 4 + base;
    }
    return index;
}

static char *read_file_contents(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) return NULL;

    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    fseek(file, 0, SEEK_SET);

    char *buffer = malloc(size + 1);
    if (!buffer) {
        fclose(file);
        return NULL;
    }

    fread(buffer, 1, size, file);
    buffer[size] = '\0';
    fclose(file);
    return buffer;
}

/* --------------- JSON Model Parsing --------------- */

int parse_json_model(const char *filename, Lambda *l, Explicit_duration *ed) {
    if (DEBUG) printf("Parsing JSON model: %s\n", filename);

    char *json_str = read_file_contents(filename);
    if (!json_str) {
        fprintf(stderr, "ERROR: Cannot read model file %s\n", filename);
        return -1;
    }

    cJSON *root = cJSON_Parse(json_str);
    if (!root) {
        fprintf(stderr, "ERROR: JSON parse error: %s\n", cJSON_GetErrorPtr());
        free(json_str);
        return -1;
    }

    /* --------------- Parse Donor PWM --------------- */
    // JSON values are already in log-space
    cJSON *don = cJSON_GetObjectItem(root, "don");
    if (don && cJSON_IsArray(don)) {
        l->B.don_kmer_len = cJSON_GetArraySize(don);
        l->B.dons         = malloc(l->B.don_kmer_len * sizeof(double*));

        for (int i = 0; i < l->B.don_kmer_len; i++) {
            l->B.dons[i]    = calloc(4, sizeof(double));
            cJSON *pos      = cJSON_GetArrayItem(don, i);

            l->B.dons[i][0] = cJSON_GetObjectItem(pos, "A")->valuedouble;
            l->B.dons[i][1] = cJSON_GetObjectItem(pos, "C")->valuedouble;
            l->B.dons[i][2] = cJSON_GetObjectItem(pos, "G")->valuedouble;
            l->B.dons[i][3] = cJSON_GetObjectItem(pos, "T")->valuedouble;
        }
        if (DEBUG) printf("  Donor PWM:    %d positions\n", l->B.don_kmer_len);
    }

    /* --------------- Parse Acceptor PWM --------------- */
    cJSON *acc = cJSON_GetObjectItem(root, "acc");
    if (acc && cJSON_IsArray(acc)) {
        l->B.acc_kmer_len = cJSON_GetArraySize(acc);
        l->B.accs         = malloc(l->B.acc_kmer_len * sizeof(double*));

        for (int i = 0; i < l->B.acc_kmer_len; i++) {
            l->B.accs[i]    = calloc(4, sizeof(double));
            cJSON *pos      = cJSON_GetArrayItem(acc, i);

            l->B.accs[i][0] = cJSON_GetObjectItem(pos, "A")->valuedouble;
            l->B.accs[i][1] = cJSON_GetObjectItem(pos, "C")->valuedouble;
            l->B.accs[i][2] = cJSON_GetObjectItem(pos, "G")->valuedouble;
            l->B.accs[i][3] = cJSON_GetObjectItem(pos, "T")->valuedouble;
        }
        if (DEBUG) printf("  Acceptor PWM: %d positions\n", l->B.acc_kmer_len);
    }

    /* --------------- Parse Exon Markov Model --------------- */
    cJSON *exs = cJSON_GetObjectItem(root, "exs");
    if (exs) {
        cJSON *k  = cJSON_GetObjectItem(exs, "k");
        cJSON *mm = cJSON_GetObjectItem(exs, "mm");

        if (k && mm) {
            l->B.exon_kmer_len = k->valueint;
            int size           = power(4, l->B.exon_kmer_len);
            l->B.exon          = calloc(size, sizeof(double));

            cJSON *entry = mm->child;
            while (entry) {
                int idx = kmer_to_index(entry->string, l->B.exon_kmer_len);
                if (idx >= 0) {
                    l->B.exon[idx] = entry->valuedouble;
                }
                entry = entry->next;
            }
            if (DEBUG) printf("  Exon MM:      k=%d (%d entries)\n",
                              l->B.exon_kmer_len, size);
        }
    }

    /* --------------- Parse Intron Markov Model --------------- */
    cJSON *ins = cJSON_GetObjectItem(root, "ins");
    if (ins) {
        cJSON *k  = cJSON_GetObjectItem(ins, "k");
        cJSON *mm = cJSON_GetObjectItem(ins, "mm");

        if (k && mm) {
            l->B.intron_kmer_len = k->valueint;
            int size              = power(4, l->B.intron_kmer_len);
            l->B.intron           = calloc(size, sizeof(double));

            cJSON *entry = mm->child;
            while (entry) {
                int idx = kmer_to_index(entry->string, l->B.intron_kmer_len);
                if (idx >= 0) {
                    l->B.intron[idx] = entry->valuedouble;
                }
                entry = entry->next;
            }
            if (DEBUG) printf("  Intron MM:    k=%d (%d entries)\n",
                              l->B.intron_kmer_len, size);
        }
    }

    /* --------------- Parse Exon Length Distribution --------------- */
    cJSON *exl = cJSON_GetObjectItem(root, "exl");
    if (exl) {
        cJSON *size_obj = cJSON_GetObjectItem(exl, "size");
        cJSON *val      = cJSON_GetObjectItem(exl, "val");
        cJSON *tail     = cJSON_GetObjectItem(exl, "tail");

        if (size_obj && val) {
            ed->exon_len     = size_obj->valueint;
            ed->max_len_exon = ed->exon_len;
            ed->exon         = calloc(ed->exon_len, sizeof(double));
            ed->min_len_exon = -1;

            double tail_val = tail ? tail->valuedouble : 0.0;
            (void)tail_val;  // May be used for extrapolation

            int arr_size = cJSON_GetArraySize(val);
            for (int i = 0; i < arr_size && i < ed->exon_len; i++) {
                cJSON *item  = cJSON_GetArrayItem(val, i);
                ed->exon[i]  = item->valuedouble;

                if (ed->exon[i] > -99.0) {  // Not log-zero (-100 sentinel)
                    if (ed->min_len_exon == -1) ed->min_len_exon = i;
                    ed->max_len_exon = i;
                }
            }
            if (DEBUG) printf("  Exon length:  min=%d, max=%d\n",
                              ed->min_len_exon, ed->max_len_exon);
        }
    }

    /* --------------- Parse Intron Length Distribution --------------- */
    cJSON *inl = cJSON_GetObjectItem(root, "inl");
    if (inl) {
        cJSON *size_obj = cJSON_GetObjectItem(inl, "size");
        cJSON *val      = cJSON_GetObjectItem(inl, "val");
        cJSON *tail     = cJSON_GetObjectItem(inl, "tail");

        if (size_obj && val) {
            ed->intron_len     = size_obj->valueint;
            ed->max_len_intron = ed->intron_len;
            ed->intron         = calloc(ed->intron_len, sizeof(double));
            ed->min_len_intron = -1;

            double tail_val = tail ? tail->valuedouble : 0.0;
            (void)tail_val;  // May be used for extrapolation

            int arr_size = cJSON_GetArraySize(val);
            for (int i = 0; i < arr_size && i < ed->intron_len; i++) {
                cJSON *item    = cJSON_GetArrayItem(val, i);
                ed->intron[i]  = item->valuedouble;

                if (ed->intron[i] > -99.0) {  // Not log-zero (-100 sentinel)
                    if (ed->min_len_intron == -1) ed->min_len_intron = i;
                    ed->max_len_intron = i;
                }
            }
            if (DEBUG) printf("  Intron length: min=%d, max=%d\n",
                              ed->min_len_intron, ed->max_len_intron);
        }
    }

    cJSON_Delete(root);
    free(json_str);

    if (DEBUG) printf("JSON model parsing complete.\n");
    return 0;
}
