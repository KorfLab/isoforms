#ifndef HMM_ALGORITHM_H
#define HMM_ALGORITHM_H

#include "types.h"
#include "model_structs.h"

/* --------------- Forward Algorithm Structures --------------- */
typedef struct
{
    double **a;                         // alpha for forward algorithm
    double **basis;                     // each previous layer of calculation
    int    first_dons;                  // where the first donor site appear
} Forward_algorithm;

/* --------------- Backward Algorithm Structures --------------- */
typedef struct
{
    double **b;                         // beta for backward algorithm
    double **basis;                     // times of transition prob and emission prob
    int    last_accs;                   // where the first acceptor site appear
} Backward_algorithm;

/* --------------- Posterior Probability --------------- */
typedef struct {
    double **xi;
    double *dons_val;
    double *accs_val;
    int    *dons_bps;
    int    *accs_bps;
    int    dons;                        // number of dons
    int    accs;                        // number of accs
} Pos_prob;

/* --------------- Viterbi Algorithm --------------- */
typedef struct {
    double **v;                         // for dp of viterbi
    double exon;
    double intron;
} Vitbi_algo;

/* --------------- Forward Algorithm Functions --------------- */
void allocate_fw(Observed_events *info, Forward_algorithm *alpha, Explicit_duration *ed);
void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info);
void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed);
void free_alpha(Observed_events *info, Forward_algorithm *alpha);

/* --------------- Backward Algorithm Functions --------------- */
void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info);
void basis_bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed);
void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed);
void free_beta(Observed_events *info, Backward_algorithm *beta);

/* --------------- Posterior Probability Functions --------------- */
void allocate_pos(Pos_prob *pos, Observed_events *info);
void free_pos(Pos_prob *pos, Observed_events *info);
void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Pos_prob *pos);
void parse_splice_sites(Pos_prob *pos, Observed_events *info);
void free_splice_sites(Pos_prob *pos);

/* --------------- Viterbi Algorithm Functions --------------- */
void allocate_vit(Vitbi_algo *vit, Observed_events *info);
void free_vit(Vitbi_algo *vit, Observed_events *info);

#endif
