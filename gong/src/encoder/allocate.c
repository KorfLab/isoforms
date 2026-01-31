#include <stdio.h>
#include <stdlib.h>
#include "model.h"

/* --------------- Forward Algorithm Memory --------------- */

void allocate_fw(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed) {
    if(DEBUG)  printf("Start allocate memory for the forward algorithm:");

    // alpha->a[t][i]
    int size_array = info->T;
    alpha->a = malloc ( (size_array) * sizeof(double*) );           // [t]: specific time among all observed events

    for( int i = 0 ; i < size_array; i++ ) {                        // [0] for exon ; [1] for intron
        alpha->a[i] = malloc( HS * sizeof(double) );                // each spot is storing a(t)(m, 1)
        for (int j = 0; j < HS; j++) alpha->a[i][j] = LOG_ZERO;     // Initialize to log(0)
    }
    // alpha->basis[i][d]
    alpha->basis    = malloc( HS * sizeof(double*) );               // [0] for exon ; [1] for intron
    // max duration for exon or intron
    alpha->basis[0] = malloc( ed->max_len_exon * sizeof(double) );
    alpha->basis[1] = malloc( ed->max_len_intron * sizeof(double));

    // Initialize basis arrays to LOG_ZERO
    for (int i = 0; i < ed->max_len_exon; i++) alpha->basis[0][i] = LOG_ZERO;
    for (int i = 0; i < ed->max_len_intron; i++) alpha->basis[1][i] = LOG_ZERO;

    if(DEBUG)  printf("\tFinished\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha) {
    if(DEBUG)  printf("Clearing up forward algorithm memory:\n");

    int size_array = info->T;
    for (int i = 0; i < size_array; i++)
        free(alpha->a[i]);

    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    if(DEBUG)  printf("\tFinished\n");
}

/* --------------- Backward Algorithm Memory --------------- */

void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info) {
    if(DEBUG)  printf("Start allocate memory for the backward algorithm:");

    // allocate basis
    beta->basis    = malloc( HS * sizeof(double*) );                    // [0] for exon ; [1] for intron
    beta->basis[0] = malloc( ed->max_len_exon   * sizeof(double) );     // max duration for exon or intron
    beta->basis[1] = malloc( ed->max_len_intron * sizeof(double) );

    // Initialize basis arrays to LOG_ZERO
    for (int i = 0; i < ed->max_len_exon; i++) beta->basis[0][i] = LOG_ZERO;
    for (int i = 0; i < ed->max_len_intron; i++) beta->basis[1][i] = LOG_ZERO;

    // allocate storage array
    int size_array = info->T;
    beta->b = malloc( (size_array) * sizeof(double*) );                 // [0] for exon ; [1] for intron
    for( int i = 0 ; i < size_array; i++ ) {
        beta->b[i] = malloc( HS * sizeof(double) );
        for (int j = 0; j < HS; j++) beta->b[i][j] = LOG_ZERO;          // Initialize to log(0)
    }

    if(DEBUG)  printf("\tFinished\n");
}

void free_beta(Observed_events *info, Backward_algorithm *beta) {
    if(DEBUG)  printf("Clearing up backward algorithm memory:\n");

    int size_array = info->T;
    for (int i = 0; i < size_array; i++)
        free(beta->b[i]);

    free(beta->b);
    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    if(DEBUG)  printf("\tFinished\n");
}

/* --------------- Posterior Probability Memory --------------- */

void allocate_pos(Pos_prob *pos, Observed_events *info){
    if(DEBUG)  printf("Start Initialize Data Structure for Posterior Probability");

    int sarray  = info->T;
    pos->xi     = malloc (sarray * sizeof(double*));
    for (int i = 0 ; i < sarray; i++) {
        pos->xi[i] = malloc( HS * sizeof(double) );
        for (int j = 0; j < HS; j++) pos->xi[i][j] = LOG_ZERO;  // Initialize to log(0)
    }

    if(DEBUG)  printf("\tFinished\n");
}

void free_pos(Pos_prob *pos, Observed_events *info) {
    if(DEBUG)  printf("Clearing up Posterior Probabilities:\n");

    int T = info->T;
    for (int t = 0; t < T; t++)
        free(pos->xi[t]);
    free(pos->xi);

    if(DEBUG)  printf("\tFinished\n");
}
