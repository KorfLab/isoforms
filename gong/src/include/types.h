#ifndef HMM_TYPES_H
#define HMM_TYPES_H

/* --------------- Core Input Structure --------------- */
typedef struct
{
    char *original_sequence;            // org seq
    int T;                              // seq len
    int *numerical_sequence;            // org seq to num seq
    int flank;                          // flank size if provided
} Observed_events;

#endif
