#include <stdio.h>
#include <stdlib.h>
#include "model.h"

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
