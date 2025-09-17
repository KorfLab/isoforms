#ifndef HMM_MODEL
#define HMM_MODEL

#define HS 2                            // 1 (exon) + 1 (intron) ; 5 (donor site) + 6(acceptor site) degraded
#define DEFAULT_FLANK 99                // default flank size if not specified

extern int  DEBUG;
extern int  use_random_forest;
extern int  n_isoforms;
extern int  use_path_restriction;

/* ---------------------------------------------------- */
/* ----------- Computation Data Structure ------------- */
/* ---------------------------------------------------- */
typedef struct
{
    char *original_sequence;            // org seq
    int T;                              // seq len
    int *numerical_sequence;            // org seq to num seq
    int flank;                          // flank size if provided
} Observed_events;

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

typedef struct {
    double *dons;                       // [dons_kmer_len]
    double *accs;                       // [accs_kmer_len]
    int don_size;                       // size of dons arr
    int acc_size;                       // size of accs arr
    // for computing all possible transition prob
    double *prob;
    int    *pos;
} Transition_matrix;

typedef struct {
    Transition_matrix A;
    Emission_matrix B;
    double *pi;
    double *log_values;
    int log_values_len;
} Lambda;

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

typedef struct
{
    double **a;                         // alpha for forward algorithm
    double **basis;                     // each previous layer of calculation
    int    first_dons;                  // where the first donor site appear
} Forward_algorithm;

typedef struct
{
    double **b;                         // beta for backward algorithm
    double **basis;                     // times of transition prob and emission prob
    int    last_accs;                   // where the first acceptor site appear
} Backward_algorithm;                   

typedef struct {
    double **xi;
    double *dons_val;
    double *accs_val;
    int    *dons_bps;
    int    *accs_bps;
    int    dons;                        // number of dons
    int    accs;                        // number of accs
} Pos_prob;

typedef struct {
    double **v;                         // for dp of viterbi
    double exon;
    double intron;
} Vitbi_algo;

/* ---------------------------------------------------- */
/* --------------- Function Declaration --------------- */
/* ---------------------------------------------------- */

/* --------------- Input Parser --------------- */
void read_sequence_file(const char *filename, Observed_events *info);
void numerical_transcription(Observed_events *info, const char *seq);
void donor_parser(Lambda *l, char *filename);
void acceptor_parser(Lambda *l, char *filename);
void exon_intron_parser(Lambda *l, char *filename, int digit);
void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit);
void compute_transition_matrices(Lambda *l);

void free_lambda(Lambda *l);
void free_explicit_duration(Explicit_duration *ed);

/* --------------- Computation Functions --------------- */
int    power(int base, int exp);
int    base4_to_int(int *array, int beg, int length);
double total_prob(double *array, int length);
double log_sum_exp(double *logs, int n);
void   tolerance_checker(double *array, int len, const double epsilon);
void   log_space_converter(double *array, int len);

/* --------------- Forward Algorithm --------------- */
void allocate_fw(Observed_events *info, Forward_algorithm *alpha, Explicit_duration *ed);
void basis_fw_algo(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info);
void fw_algo(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed);
void free_alpha(Observed_events *info, Forward_algorithm *alpha);

/* --------------- Backward Algorithm --------------- */
void allocate_bw(Backward_algorithm *beta, Explicit_duration *ed, Observed_events *info);
void basis_bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed);
void bw_algo(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed);
void free_beta(Observed_events *info, Backward_algorithm *beta);

/* --------------- Posterior Probability --------------- */
void allocate_pos(Pos_prob *pos, Observed_events *info);
void free_pos(Pos_prob *pos, Observed_events *info);
void pos_prob(Backward_algorithm *beta, Forward_algorithm *alpha, Observed_events *info, Pos_prob *pos);
void parse_splice_sites(Pos_prob *pos, Observed_events *info);
void free_splice_sites(Pos_prob *pos);

/* --------------- Output --------------- */
void index_to_sequence(int index, int length, char *seq);
void print_transition_matrices_summary(Lambda *l);
void print_splice_sites(Pos_prob *pos, Observed_events *info);
void print_duration_summary(Explicit_duration *ed);

/* --------------- Viterbi Algorithm Functions --------------- */
void allocate_vit(Vitbi_algo *vit, Observed_events *info);
void free_vit(Vitbi_algo *vit, Observed_events *info);

#endif