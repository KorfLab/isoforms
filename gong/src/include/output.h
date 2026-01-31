#ifndef HMM_OUTPUT_H
#define HMM_OUTPUT_H

#include "types.h"
#include "model_structs.h"
#include "algorithm.h"

/* --------------- Debug Output --------------- */
void index_to_sequence(int index, int length, char *seq);
void print_transition_matrices_summary(Lambda *l);
void print_duration_summary(Explicit_duration *ed);

/* --------------- Splice Site Output --------------- */
void print_splice_sites(Pos_prob *pos, Observed_events *info);

#endif
