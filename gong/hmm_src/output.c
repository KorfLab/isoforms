#include "model.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ==================================================== *
 * =============== HMM Hints_output =================== *
 * ==================================================== */

void print_splice_sites(Pos_prob *pos, Observed_events *info, Explicit_duration *ed)
{
    printf("DONS\n");
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][0] != 0.0) {
            printf("%d\t%.10f\n", i, pos->xi[i][0]);
        }
    }
    
    printf("ACCS\n");
    for (int i = FLANK; i < info->T-FLANK; i++) {
        if (pos->xi[i][1] != 0.0) {  
            printf("%d\t%.10f\n", i, pos->xi[i][1]);
        }
    }
}

/* ==================================================== *
 * ============== Parser Debug Output ================= *
 * ==================================================== */

void index_to_sequence(int index, int length, char *seq)
{
    for(int i = length-1; i >= 0; i--) {
        int base = index % 4;
        switch(base) {
            case 0: seq[i] = 'A'; break;
            case 1: seq[i] = 'C'; break;
            case 2: seq[i] = 'G'; break;
            case 3: seq[i] = 'T'; break;
        }
        index /= 4;
    }
    seq[length] = '\0';
}

void print_transition_matrices_summary(Lambda *l)
{
    printf("\n=== TRANSITION MATRIX SUMMARY ===\n");
    
    // Donor summary
    int non_zero_donors = 0;
    double max_donor = 0.0;
    int max_donor_idx = -1;
    
    for(int i = 0; i < 1024; i++) {
        double prob = l->A.dons[i];
        if(prob > 0) {
            non_zero_donors++;
            if(prob > max_donor) {
                max_donor = prob;
                max_donor_idx = i;
            }
        }
    }
    
    printf("Donor Summary:\n");
    printf("Non-zero entries: %d / 1024\n", non_zero_donors);
    if(max_donor_idx >= 0) {
        char max_seq[6];
        index_to_sequence(max_donor_idx, 5, max_seq);
        printf("Maximum probability: %.6e (sequence: %s, index: %d)\n", max_donor, max_seq, max_donor_idx);
    }
    
    // Acceptor summary  
    int non_zero_acceptors = 0;
    double max_acceptor = 0.0;
    int max_acceptor_idx = -1;
    
    for(int i = 0; i < 4096; i++) {
        double prob = l->A.accs[i];
        if(prob > 0) {
            non_zero_acceptors++;
            if(prob > max_acceptor) {
                max_acceptor = prob;
                max_acceptor_idx = i;
            }
        }
    }
    
    printf("Acceptor Summary:\n");
    printf("Non-zero entries: %d / 4096\n", non_zero_acceptors);
    if(max_acceptor_idx >= 0) {
        char max_seq[7];
        index_to_sequence(max_acceptor_idx, 6, max_seq);
        printf("Maximum probability: %.6e (sequence: %s, index: %d)\n", max_acceptor, max_seq, max_acceptor_idx);
    }
    
    printf("================================\n\n");
}

void print_duration_summary(Explicit_duration *ed)
{
    printf("\n=== EXPLICIT DURATION SUMMARY ===\n");
    
    // Exon analysis
    printf("EXON DURATION:\n");
    printf("  Min length: %d, Max length: %d\n", ed->min_len_exon, ed->max_len_exon);
    
    int exon_nonzero = 0;
    int exon_first_nonzero = -1;
    int exon_last_nonzero = -1;
    double exon_sum = 0.0;
    double exon_max_prob = 0.0;
    int exon_max_idx = -1;
    
    for(int i = 0; i < ed->max_len_exon; i++) {
        if(ed->exon[i] > 0.0) {
            exon_nonzero++;
            exon_sum += ed->exon[i];
            
            if(exon_first_nonzero == -1) exon_first_nonzero = i;
            exon_last_nonzero = i;
            
            if(ed->exon[i] > exon_max_prob) {
                exon_max_prob = ed->exon[i];
                exon_max_idx = i;
            }
        }
    }
    
    printf("  Non-zero entries: %d / %d\n", exon_nonzero, ed->max_len_exon);
    printf("  First non-zero at: %d (should match min_len: %d) %s\n", 
           exon_first_nonzero, ed->min_len_exon, 
           (exon_first_nonzero == ed->min_len_exon) ? "✓" : "✗");
    printf("  Last non-zero at: %d\n", exon_last_nonzero);
    printf("  Range span: %d (last - first + 1 = %d)\n", 
           exon_last_nonzero - exon_first_nonzero + 1, exon_nonzero);
    printf("  Probability sum: %.6f (should be ~1.0) %s\n", 
           exon_sum, (fabs(exon_sum - 1.0) < 0.01) ? "✓" : "✗");
    printf("  Maximum prob: %.6e at length %d\n", exon_max_prob, exon_max_idx);
    
    // Intron analysis
    printf("\nINTRON DURATION:\n");
    printf("  Min length: %d, Max length: %d\n", ed->min_len_intron, ed->max_len_intron);
    
    int intron_nonzero = 0;
    int intron_first_nonzero = -1;
    int intron_last_nonzero = -1;
    double intron_sum = 0.0;
    double intron_max_prob = 0.0;
    int intron_max_idx = -1;
    
    for(int i = 0; i < ed->max_len_intron; i++) {
        if(ed->intron[i] > 0.0) {
            intron_nonzero++;
            intron_sum += ed->intron[i];
            
            if(intron_first_nonzero == -1) intron_first_nonzero = i;
            intron_last_nonzero = i;
            
            if(ed->intron[i] > intron_max_prob) {
                intron_max_prob = ed->intron[i];
                intron_max_idx = i;
            }
        }
    }
    
    printf("  Non-zero entries: %d / %d\n", intron_nonzero, ed->max_len_intron);
    printf("  First non-zero at: %d (should match min_len: %d) %s\n", 
           intron_first_nonzero, ed->min_len_intron, 
           (intron_first_nonzero == ed->min_len_intron) ? "✓" : "✗");
    printf("  Last non-zero at: %d\n", intron_last_nonzero);
    printf("  Range span: %d (last - first + 1 = %d)\n", 
           intron_last_nonzero - intron_first_nonzero + 1, intron_nonzero);
    printf("  Probability sum: %.6f (should be ~1.0) %s\n", 
           intron_sum, (fabs(intron_sum - 1.0) < 0.01) ? "✓" : "✗");
    printf("  Maximum prob: %.6e at length %d\n", intron_max_prob, intron_max_idx);
    
    // Overall validation
    printf("\nVALIDATION:\n");
    int errors = 0;
    
    if(exon_first_nonzero != ed->min_len_exon) {
        printf("  ✗ Exon min_len mismatch!\n");
        errors++;
    }
    if(intron_first_nonzero != ed->min_len_intron) {
        printf("  ✗ Intron min_len mismatch!\n");
        errors++;
    }
    if(fabs(exon_sum - 1.0) > 0.01) {
        printf("  ✗ Exon probabilities don't sum to 1.0!\n");
        errors++;
    }
    if(fabs(intron_sum - 1.0) > 0.01) {
        printf("  ✗ Intron probabilities don't sum to 1.0!\n");
        errors++;
    }
    
    if(errors == 0) {
        printf("  ✓ All duration probabilities look good!\n");
    } else {
        printf("  Found %d validation errors!\n", errors);
    }
    
    printf("==================================\n\n");
}