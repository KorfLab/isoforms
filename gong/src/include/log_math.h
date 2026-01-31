#ifndef HMM_LOG_MATH_H
#define HMM_LOG_MATH_H

#include <math.h>
#include "constants.h"

/* --------------- Log-Space Helper Macros --------------- */
// Check if a log-space value represents zero probability
#define IS_LOG_ZERO(x) ((x) <= LOG_ZERO + 1.0)

// Safe log: returns LOG_ZERO for non-positive values
#define SAFE_LOG(x) ((x) > 0.0 ? log(x) : LOG_ZERO)

// Safe exp: prevents overflow/underflow
#define SAFE_EXP(x) ((x) > LOG_EPSILON ? exp(x) : 0.0)

// Log-space multiplication: simply addition in log space
#define LOG_MUL(a, b) ((IS_LOG_ZERO(a) || IS_LOG_ZERO(b)) ? LOG_ZERO : ((a) + (b)))

/* --------------- Log-Space Inline Functions --------------- */
// Log-space addition: log(exp(a) + exp(b)) with numerical stability
static inline double log_add(double a, double b) {
    if (IS_LOG_ZERO(a)) return b;
    if (IS_LOG_ZERO(b)) return a;
    if (a > b) return a + log1p(exp(b - a));
    else return b + log1p(exp(a - b));
}

/* --------------- Utility Functions --------------- */
int    power(int base, int exp);
int    base4_to_int(int *array, int beg, int length);
double total_prob(double *array, int length);
double log_sum_exp(double *logs, int n);
void   tolerance_checker(double *array, int len, const double epsilon);
void   log_space_converter(double *array, int len);

#endif
