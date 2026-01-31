#ifndef HMM_CONSTANTS_H
#define HMM_CONSTANTS_H

/* --------------- Model Constants --------------- */
#define HS 2                            // 1 (exon) + 1 (intron) ; 5 (donor site) + 6(acceptor site) degraded
#define DEFAULT_FLANK 99                // default flank size if not specified

/* --------------- Numerical Stability Constants --------------- */
#define LOG_ZERO (-1e300)               // Represents log(0) = -infinity, but finite for computation
#define LOG_ONE 0.0                     // log(1) = 0
#define LOG_EPSILON (-700.0)            // Minimum log value before underflow (exp(-700) ≈ 0)
#define PROB_EPSILON 1e-300             // Minimum probability before treating as zero

/* --------------- Global Flags --------------- */
extern int  DEBUG;
extern int  use_random_forest;
extern int  n_isoforms;
extern int  use_path_restriction;
extern int  model_in_log_space;

#endif
