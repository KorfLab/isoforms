#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "model.h"
#include "decoder/randomf.h"

int     DEBUG                   = 0;
int     use_random_forest       = 0;
int     n_isoforms              = 10000;
int     model_in_log_space      = 0;

void print_usage(const char *program_name) {
    printf("RFHMM - Random Forest Hidden Markov Model for gene prediction\n");
    printf("Usage: %s [OPTIONS]\n\n", program_name);
    printf("Required options:\n");
    printf("  -s, --sequence FILE           Input sequence file\n");
    printf("\nModel options (use ONE of the following):\n");
    printf("  -M, --model FILE              JSON model file (.splicemodel)\n");
    printf("  OR individual model files:\n");
    printf("  -d, --don_emission FILE       Donor emission file (default: ../models/don.pwm)\n");
    printf("  -a, --acc_emission FILE       Acceptor emission file (default: ../models/acc.pwm)\n");
    printf("  -e, --exon_emission FILE      Exon emission file (default: ../models/exon.mm)\n");
    printf("  -i, --intron_emission FILE    Intron emission file (default: ../models/intron.mm)\n");
    printf("  -x, --ped_exon FILE           Exon length distribution file (default: ../models/exon.len)\n");
    printf("  -n, --ped_intron FILE         Intron length distribution file (default: ../models/intron.len)\n");
    printf("\nAnalysis parameters:\n");
    printf("  -f, --flank NUM               Flank size for analysis (default: 99)\n");
    printf("\nRandom Forest options:\n");
    printf("  -S, --stovit                  Enable stochastic Viterbi with Random Forest\n");
    printf("  -N, --n_isoforms NUM          Maximum isoform capacity (default: 10000)\n");
    printf("  -m, --mtry FRACTION           Fraction of features to sample (e.g., 0.5 or 1/2, default: 0.5)\n");
    printf("  -z, --node_size NUM           Minimum number of observations in terminal node (default: 5)\n");
    printf("\nOutput control:\n");
    printf("  -p, --print_splice            Print detailed splice site analysis\n");
    printf("  -j, --json                    Emit isoform locus information as JSON\n");
    printf("  -v, --verbose                 Show debug and progress information\n");
    printf("  -h, --help                    Show this help message\n");
    printf("\nExamples:\n");
    printf("  %s --sequence input.fasta \n", program_name);
    printf("  %s -s input.fasta --stovit\n", program_name);
    printf("  %s -s input.fasta --stovit --n_isoforms 5000 --mtry 3\n", program_name);
    printf("  %s -s input.fasta --print_splice --flank 150\n", program_name);
    printf("  %s -s input.fasta --verbose\n", program_name);
    printf("\nRandom Forest Algorithm:\n");
    printf("  The algorithm generates trees continuously until the locus capacity\n");
    printf("  is reached. It automatically finds unique isoforms without needing\n");
    printf("  to specify target counts or tree limits.\n");
}

int main(int argc, char *argv[])
{
    // Hard Code Path for Default(shall be deleted)
    char *default_don_emission      = "../models/don.pwm";
    char *default_acc_emission      = "../models/acc.pwm";
    char *default_exon_emission     = "../models/exon.mm";
    char *default_intron_emission   = "../models/intron.mm";
    char *default_Ped_exon          = "../models/exon.len";
    char *default_Ped_intron        = "../models/intron.len";

    // Variables for command-line inputs
    char *don_emission          = default_don_emission;
    char *acc_emission          = default_acc_emission;
    char *exon_emission         = default_exon_emission;
    char *intron_emission       = default_intron_emission;
    char *Ped_exon              = default_Ped_exon;
    char *Ped_intron            = default_Ped_intron;
    char *seq_input             = NULL;
    char *model_file            = NULL;
    int print_splice_detailed   = 0;
    int output_json             = 0;
    int flank_size              = DEFAULT_FLANK;
    float mtry                  = 0.5;          // Default to 1/2 of features
    int node_size               = 5;            // for regression 5 as node_size

    static struct option long_options[] = {
        {"sequence",        required_argument, 0, 's'},
        {"model",           required_argument, 0, 'M'},
        {"don_emission",    required_argument, 0, 'd'},
        {"acc_emission",    required_argument, 0, 'a'},
        {"exon_emission",   required_argument, 0, 'e'},
        {"intron_emission", required_argument, 0, 'i'},
        {"ped_exon",        required_argument, 0, 'x'},
        {"ped_intron",      required_argument, 0, 'n'},
        {"flank",           required_argument, 0, 'f'},
        {"n_isoforms",      required_argument, 0, 'N'},
        {"mtry",            required_argument, 0, 'm'},
        {"node_size",       required_argument, 0, 'z'},
        {"json",            no_argument,       0, 'j'},
        {"print_splice",    no_argument,       0, 'p'},
        {"stovit",          no_argument,       0, 'S'},
        {"verbose",         no_argument,       0, 'v'},
        {"help",            no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "s:M:d:a:e:i:x:n:f:N:m:z:jSpvh", long_options, &option_index)) != -1) {
        switch (c) {
            case 's':
                seq_input = optarg;
                break;
            case 'M':
                model_file = optarg;
                break;
            case 'd':
                don_emission = optarg;
                break;
            case 'a':
                acc_emission = optarg;
                break;
            case 'e':
                exon_emission = optarg;
                break;
            case 'i':
                intron_emission = optarg;
                break;
            case 'x':
                Ped_exon = optarg;
                break;
            case 'n':
                Ped_intron = optarg;
                break;
            case 'f':
                flank_size = atoi(optarg);
                if (flank_size < 0) {
                    fprintf(stderr, "Error: flank size must be non-negative\n");
                    return 1;
                }
                break;
            case 'p':
                print_splice_detailed = 1;
                break;
            case 'v':
                DEBUG = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'N':
                n_isoforms = atoi(optarg);
                if (n_isoforms < 1) {
                    fprintf(stderr, "Error: n_isoforms must be at least 1\n");
                    return 1;
                }
                break;
            case 'm':
                if (strchr(optarg, '/')) {
                    int numerator, denominator;
                    if (sscanf(optarg, "%d/%d", &numerator, &denominator) == 2 && denominator != 0) {
                        mtry = (float)numerator / denominator;
                    } else {
                        fprintf(stderr, "Error: Invalid fraction format for mtry\n");
                        return 1;
                    }
                } else {
                    mtry = atof(optarg);
                }
                if (mtry <= 0.0 || mtry > 1.0) {
                    fprintf(stderr, "Error: mtry must be between 0 and 1 (got %.2f)\n", mtry);
                    return 1;
                }
                break;
            case 'z':
                node_size = atoi(optarg);
                if (node_size < 1) {
                    fprintf(stderr, "Error: node_size must be at least 1\n");
                    return 1;
                }
                break;
            case 'j':
                output_json = 1;
                break;
            case 'S':
                use_random_forest = 1;
                break;
            case '?':
                // getopt_long already printed an error message
                print_usage(argv[0]);
                return 1;
            default:
                abort();
        }
    }

    if (seq_input == NULL) {
        fprintf(stderr, "Error: --sequence is required\n");
        print_usage(argv[0]);
        return 1;
    }

    /* --------------- Initialize Data Structure --------------- */
    Observed_events info;
    Lambda l;
    Explicit_duration ed;
    Forward_algorithm fw;
    Backward_algorithm bw;
    Pos_prob pos;

    memset(&info, 0, sizeof(Observed_events));
    info.flank = flank_size;
    memset(&ed, 0, sizeof(Explicit_duration));
    ed.min_len_exon     = -1;
    ed.min_len_intron   = -1;
    ed.max_len_exon     = 0;
    ed.max_len_intron   = 0;
    memset(&l, 0, sizeof(Lambda));

    if (DEBUG) printf("\n=== Starting RFHMM Analysis ===\n");

    /* --------------- Parse Sequence --------------- */
    if (DEBUG) printf("\n--- Phase 1: Loading Sequence ---\n");
    read_sequence_file(seq_input, &info);
    numerical_transcription(&info, info.original_sequence);
    
    // FLANK Size Check
    if (info.flank * 2 >= info.T) {
        fprintf(stderr, "Error: Flank size (%d) is too large for sequence length (%d)\n", 
                info.flank, info.T);
        fprintf(stderr, "Flank size should be less than %d\n", info.T / 2);
        return 1;
    }

    /* --------------- Parse HMM Model Input --------------- */
    if (DEBUG) printf("\n--- Phase 2: Parsing Model Files ---\n");

    if (model_file) {
        // Use unified JSON model file (values already in log-space)
        if (parse_json_model(model_file, &l, &ed) != 0) {
            fprintf(stderr, "Error: Failed to parse model file %s\n", model_file);
            return 1;
        }
        model_in_log_space = 1;
    } else {
        // Parse individual model files
        donor_parser(&l, don_emission);
        acceptor_parser(&l, acc_emission);
        exon_intron_parser(&l, exon_emission, 0);
        exon_intron_parser(&l, intron_emission, 1);

        explicit_duration_probability(&ed, Ped_exon, 0);
        explicit_duration_probability(&ed, Ped_intron, 1);
    }

    if (DEBUG) printf("\n--- Phase 3: Computing Transition Matrices ---\n");
    if (model_in_log_space) {
        compute_transition_matrices_log(&l);
    } else {
        compute_transition_matrices(&l);
    }
    
    l.log_values_len = (ed.max_len_exon > ed.max_len_intron) ? ed.max_len_exon : ed.max_len_intron;
    l.log_values     = calloc(l.log_values_len, sizeof(double));

    /* --------------- Validation Check --------------- */
    if (DEBUG) {
        printf("\n--- Phase 4: Validation ---\n");
        printf("Parameters loaded:\n");
        printf("  Sequence length: %d bp\n", info.T);
        printf("  Flank size: %d bp\n", info.flank);
        printf("  Exon length range: %d-%d bp\n", ed.min_len_exon, ed.max_len_exon);
        printf("  Intron length range: %d-%d bp\n", ed.min_len_intron, ed.max_len_intron);
        printf("  Analysis range: %d to %d\n", 
               info.flank+ed.min_len_exon, info.T-info.flank-ed.min_len_exon);
        
        print_transition_matrices_summary(&l);
        print_duration_summary(&ed);
    }

    /* --------------- Log Space Conversion --------------- */
    int don_size    = power(4, l.B.don_kmer_len);
    int acc_size    = power(4, l.B.acc_kmer_len);
    int exon_size   = power(4, l.B.exon_kmer_len);
    int intron_size = power(4, l.B.intron_kmer_len);

    if (model_in_log_space) {
        // JSON model already in log-space, skip conversion
        if (DEBUG) printf("\n--- Phase 5: Model already in log space ---\n");
    } else {
        if (DEBUG) printf("\n--- Phase 5: Converting to Log Space ---\n");

        // Check for numerical issues
        tolerance_checker(ed.exon, ed.exon_len, 1e-15);
        tolerance_checker(ed.intron, ed.intron_len, 1e-15);
        tolerance_checker(l.A.dons, don_size, 1e-15);
        tolerance_checker(l.A.accs, acc_size, 1e-15);

        // Convert to log space for numerical stability
        log_space_converter(ed.exon, ed.exon_len);
        log_space_converter(ed.intron, ed.intron_len);
        log_space_converter(l.A.dons, don_size);
        log_space_converter(l.A.accs, acc_size);
        log_space_converter(l.B.exon, exon_size);
        log_space_converter(l.B.intron, intron_size);
    }

    /* --------------- Exe Forward Backward Algorithm --------------- */
    if (DEBUG) printf("\n--- Phase 6: Forward-Backward Algorithm ---\n");
    
    // Allocate memory for algorithms
    allocate_fw(&info, &fw, &ed);
    allocate_bw(&bw, &ed, &info);
    allocate_pos(&pos, &info);

    // Run forward algorithm
    basis_fw_algo(&l, &ed, &fw, &info);
    if (DEBUG) printf("  Forward basis complete\n");
    fw_algo(&l, &fw, &info, &ed);
    if (DEBUG) printf("  Forward algorithm complete\n");
    
    // Run backward algorithm
    basis_bw_algo(&l, &bw, &info, &ed);
    if (DEBUG) printf("  Backward basis complete\n");
    bw_algo(&l, &bw, &info, &ed);
    if (DEBUG) printf("  Backward algorithm complete\n");
    
    // Calculate posterior probabilities
    pos_prob(&bw, &fw, &info, &pos);
    if (DEBUG) printf("  Posterior probabilities calculated\n");
    
    // Parse splice sites from posterior probabilities
    parse_splice_sites(&pos, &info);
    
    if (DEBUG) {
        printf("\n=== Results Summary ===\n");
        printf("Found %d donor sites and %d acceptor sites\n", pos.dons, pos.accs);
    }

    /* --------------- For Stovit --------------- */
    if (use_random_forest) {
        if (DEBUG) printf("\n--- Phase 7: Random Forest Isoform Generation ---\n");
        
        if (pos.dons == 0 || pos.accs == 0) {
            printf("Warning: No splice sites found. Cannot run Random Forest.\n");
        } else {
            Locus *loc = create_locus(n_isoforms);
            Vitbi_algo vit;
            memset(&vit, 0, sizeof(Vitbi_algo));
            allocate_vit(&vit, &info);
            
            if (DEBUG) {
                printf("Generating isoforms using Random Forest:\n");
                printf("  Locus capacity: %d\n", n_isoforms);
                printf("  Node size: %d\n", node_size);
                printf("  Mtry: %.2f (%.0f%% of features)\n", mtry, mtry * 100);
                printf("  Path restriction: Yes\n");
            }

            RandomForest *rf = create_random_forest(&pos, loc, node_size, mtry);
            generate_isoforms_random_forest(rf, &info, &ed, &l, loc, &vit);
            if (DEBUG) printf("Unique isoforms found: %d\n", loc->n_isoforms);
            if (output_json) {
                print_locus_json(loc, &info, stdout);
            } else {
                print_locus(loc, &info);    
            }

            // Clean up random forest and locus
            free_random_forest(rf);
            free_vit(&vit, &info);
            free_locus(loc);
        }
    }

    /* --------------- For HMM Hints --------------- */
    if (print_splice_detailed) {
        print_splice_sites(&pos, &info);
    }

    /* --------------- Memory Cleanup --------------- */
    if (DEBUG) printf("\n--- Phase 8: Cleanup ---\n");
    
    free_splice_sites(&pos);
    free_alpha(&info, &fw);
    free_beta(&info, &bw);
    free_pos(&pos, &info);
    free(info.original_sequence);
    free(info.numerical_sequence);
    free_lambda(&l);
    free_explicit_duration(&ed);
    
    if (DEBUG) printf("\n=== RFHMM Analysis Complete ===\n");
    return 0;
}