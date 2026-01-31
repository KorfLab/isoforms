#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

/* --------------- Sequence Transcription --------------- */

void numerical_transcription(Observed_events *info, const char *seq) {
    if(DEBUG)  printf("Start transforming original sequence into base4:\n");
    // turns original sequence into int
    size_t len  = strlen(seq);
    info->T     = len;
    info->numerical_sequence = malloc ( len * sizeof(int) );

    for (size_t i = 0; i < len; i++) {
        if      (seq[i] == 'A')     info->numerical_sequence[i] = 0;
        else if (seq[i] == 'C')     info->numerical_sequence[i] = 1;
        else if (seq[i] == 'G')     info->numerical_sequence[i] = 2;
        else if (seq[i] == 'T')     info->numerical_sequence[i] = 3;
    }

    if(DEBUG)  printf("\tWe get numerical sequence with Seq len: %d\n", info->T);
    if(DEBUG)  printf("\tFinished\n");
    if(DEBUG)  printf("\n");
}

/* --------------- Math Utilities --------------- */

int power(int base, int exp) {
    int result = 1;
    for (int i = 0 ; i < exp ; i++) {
        result *= base;
    }
    return result;
}

int base4_to_int(int *array, int beg, int length) {
    int values = 0;
    int value;
    for (int i = 0; i < length; i++) {
        value  =  array[beg + i];
        values += value * power(4, length - i - 1);
    }
    return values;
}

double log_sum_exp(double *array, int n) {
    if(!array){
        printf("Something wrong with log_sum_exp trick. Invalid Input\n");
        return LOG_ZERO;
    }
    if(n <= 0) return LOG_ZERO;

    // Find max, skipping LOG_ZERO values
    double max = LOG_ZERO;
    int valid_count = 0;
    for (int i = 0; i < n; i++) {
        if (!IS_LOG_ZERO(array[i])) {
            if (valid_count == 0 || array[i] > max) {
                max = array[i];
            }
            valid_count++;
        }
    }

    if (valid_count == 0) return LOG_ZERO;

    // Compute sum in numerically stable way
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        if (!IS_LOG_ZERO(array[i])) {
            sum += exp(array[i] - max);
        }
    }
    return max + log(sum);
}

/* --------------- Validators --------------- */

void tolerance_checker(double *array, int len, const double epsilon) {
    if(!array){
        printf("This is not a valid input for float20.0\n");
        return;
    }

    for(int i = 0 ; i < len ; i++) {
        if( array[i] < epsilon )    array[i] = 0;
    }
}

void log_space_converter(double *array, int len) {
    if(!array){
        printf("This is not a valid input for log_space converter\n");
        return;
    }

    for (int i = 0; i < len ; i++) {
        if (array[i] <= PROB_EPSILON)   array[i] = LOG_ZERO;
        else                            array[i] = log(array[i]);
    }
}
