#ifndef HMM_PARSER_H
#define HMM_PARSER_H

#include "types.h"
#include "model_structs.h"

/* --------------- Sequence Parsing --------------- */
void read_sequence_file(const char *filename, Observed_events *info);
void numerical_transcription(Observed_events *info, const char *seq);

/* --------------- Model File Parsing --------------- */
void donor_parser(Lambda *l, char *filename);
void acceptor_parser(Lambda *l, char *filename);
void exon_intron_parser(Lambda *l, char *filename, int digit);
void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit);

/* --------------- Transition Matrix Computation --------------- */
void compute_transition_matrices(Lambda *l);
void compute_transition_matrices_log(Lambda *l);

/* --------------- JSON Model Parsing --------------- */
int  parse_json_model(const char *filename, Lambda *l, Explicit_duration *ed);

/* --------------- Memory Management --------------- */
void free_lambda(Lambda *l);
void free_explicit_duration(Explicit_duration *ed);

#endif
