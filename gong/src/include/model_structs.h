#ifndef HMM_MODEL_STRUCTS_H
#define HMM_MODEL_STRUCTS_H

/* --------------- Emission Matrix --------------- */
typedef struct {
    double **dons;                      // [dons kmer len][4]
    double **accs;                      // [accs kmer len][4]
    double *exon;                       // [exon   kmer len]
    double *intron;                     // [intron kmer len]
    int don_kmer_len;                   // dons     k-mer len
    int acc_kmer_len;                   // accs     k-mer len
    int exon_kmer_len;                  // exon     k-mer len
    int intron_kmer_len;                // intron   k-mer len
} Emission_matrix;

/* --------------- Transition Matrix --------------- */
typedef struct {
    double *dons;                       // [dons_kmer_len]
    double *accs;                       // [accs_kmer_len]
    int don_size;                       // size of dons arr
    int acc_size;                       // size of accs arr
    // for computing all possible transition prob
    double *prob;
    int    *pos;
} Transition_matrix;

/* --------------- HMM Parameters --------------- */
typedef struct {
    Transition_matrix A;
    Emission_matrix B;
    double *pi;
    double *log_values;
    int log_values_len;
} Lambda;

/* --------------- Explicit Duration --------------- */
typedef struct {
    double *exon;
    double *intron;
    int exon_len;                       // len of exon   arr
    int intron_len;                     // len of intron arr
    int max_len_exon;
    int max_len_intron;
    int min_len_exon;
    int min_len_intron;
} Explicit_duration;

#endif
