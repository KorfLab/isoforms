#ifndef HMM_MODEL
#define HMM_MODEL

#define HS 2                            // 1 (exon) + 1 (intron) ; 5 (donor site) + 6(acceptor site) degraded
#define FLANK 25                        // define the global flank size

typedef struct                          // observed events with length T
{
    char *original_sequence;            // where the original sequence store
    int T;                              // overall length for sequence
    int *numerical_sequence;            // transcribe from base pair to digits
} Observed_events;

typedef struct
{
    double dons[5][4];                  // the emission probability for donor sites
    double accs[6][4];                  // the emission probability for acceptor sites
    double exon[256];                   // the emission probability for exon
    double intron[256];                 // the emission probability for intron
} Emission_matrix;

typedef struct                          // degrade the sequence of conventional transition prob from donor 1-5 acceptor 1-6
{
    double dons[1024];                  // enumerating exon->intron ; aka donor site series
    double accs[4096];                  // enumerating intron->exon ; aka acceptor site series                       
} Transition_matrix;

typedef struct 
{
    double prob[6];                     // for apc algorithm to calculate transition prob
    int position[6];                    // for apc algorithm to get index to store in transition matrix
} Apc;

typedef struct
{
    Transition_matrix A;                // the transition probability
    Emission_matrix B;                  // the pre-defined emission probibility data strcuture
    double *pi;                         // the initial probability
    double log_values[999];            // prepared for log softmax trick
} Lambda;

typedef struct
{
    double exon[1000];                  // the ed probability for exon
    double intron[1000];                // the ed probability for intron
    int max_len_exon;                   // max len for exon
    int max_len_intron;                 // max len for intron
    int min_len_exon;                   // min len for exon
    int min_len_intron;                 // min len for intron
} Explicit_duration;

typedef struct
{
    double **a;                         // alpha for forward algorithm
    double **basis;                     // each previous layer of calculation
} Forward_algorithm;

typedef struct
{
    double **basis;                     // times of transition prob and emission prob
} Backward_algorithm;                   

typedef struct
{
    double *xi;
    double *gamma;
    int *path;                        
    double **xi_sum;
    double xi_sum_exon;
    double xi_sum_intron;
} Viterbi_algorithm;



// declared function //

// seq reading //

void read_sequence_file(const char *filename, Observed_events *info);
void numerical_transcription(Observed_events *info, const char *seq);

// input model //

void donor_parser(Lambda *l, char *filename);
void acceptor_parser(Lambda *l, char *filename);
void exon_intron_parser(Lambda *l, char *filename, int digit);
void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit);

// EDHMM setup // 

void setup_initial_probability(Lambda *l);

// computation function //

void normalize_transition_prob(Lambda *l, int len, int dons_or_accs);
int power(int base, int exp);
int base4_to_int(int *array, int beg, int length);
double total_prob(double *array, int length);
double safe_log(double x);
double log_sum_exp(double *logs, int n);

// suffix algorithm to get all transition prob // 
void initialize_donor_transition_matrix(Lambda *l, Apc *a, int depth);
void initialize_acceptor_transition_matrix(Lambda *l, Apc *a, int depth);

// forward algorithm //

void allocate_alpha(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed);                        
void basis_forward_algorithm(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info);
void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed);
void free_alpha(Observed_events *info, Forward_algorithm *alpha);

// viterbi algorithm //

void allocate_viterbi(Viterbi_algorithm *vit, Observed_events *info);
void viterbi_basis(Viterbi_algorithm *vit, Forward_algorithm *alpha);
void argmax_viterbi(Viterbi_algorithm *vit, int t);
void xi_calculation(Lambda *l, Forward_algorithm *alpha, Viterbi_algorithm *vit, Observed_events *info, double backward_sum, int t, int type);
void free_viterbi(Viterbi_algorithm *vit);

// backward algorithm //

void allocate_beta(Backward_algorithm *beta, Explicit_duration *ed);                             
void initial_backward_algorithm(Backward_algorithm *beta);
void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit, Forward_algorithm *alpha);
void free_beta(Backward_algorithm *beta);

// output section //
void viterbi_path_test(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed);

#endif