#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "model.h"

/* ==================================================== *
 * ==================== Seq reader ==================== *
 * ==================================================== */

void read_sequence_file(const char *filename, Observed_events *info) {
    if (DEBUG == 1 || DEBUG == 2) printf("Start reading the sequence data:\n");

    FILE *file = fopen(filename, "r");
    
    if (file == NULL) {
        if (DEBUG == 1 || DEBUG == 2) printf("Error: Cannot open sequence file %s\n", filename);
        return;
    }
    
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    char *buffer = (char*)malloc(file_size + 1);
    size_t read_size = fread(buffer, 1, file_size, file);

    buffer[read_size] = '\0';
    char *sequence = (char*)malloc(file_size + 1);

    size_t seq_index = 0;
    int in_header = 0;
    
    char *line = strtok(buffer, "\n");
    while (line != NULL) {
        if (line[0] == '>') {
            in_header = 1;
            if (DEBUG == 1 || DEBUG == 2) printf("\tSkipping header: %s\n", line);
            line = strtok(NULL, "\n");
            continue;
        }
        
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T') {
                sequence[seq_index++] = line[i];
            }
        }
        
        line = strtok(NULL, "\n");
    }
    
    sequence[seq_index] = '\0';
    sequence = (char*)realloc(sequence, seq_index + 1);
    
    info->original_sequence = sequence;
    info->T = seq_index;
    
    free(buffer);
    fclose(file);

    if (DEBUG == 1 || DEBUG == 2) printf("\tWe get original sequence with Seq len: %zu\n", seq_index);
    if (DEBUG == 1 || DEBUG == 2) printf("\tFinished\n");
    if (DEBUG == 1 || DEBUG == 2) printf("\n");
}

/* ==================================================== *
 * ============== Transition Probability ============== *
 * ==================================================== */

void donor_parser(Lambda *l, char *filename)            // get emission probability for donor site
{
    if (DEBUG == 1 || DEBUG == 2)     printf("Start getting donor site emission Probability:");
    FILE *file = fopen(filename, "r");

    char line[256];
    char *token;
    double p;                                           // probability we are going to store

    if (file == NULL)
    {
        if (DEBUG == 1 || DEBUG == 2)     printf("Can't find file for donor site emission probability!\n");
        return;
    }
    
    int c_line = -1;                                    // count of line

    while( fgets( line, sizeof(line) , file) != NULL )  // nest while loop to get elements
    {

        if ( line[0] == '%')     continue;              // skip the first line        

        c_line++;
        int c_token = -1;

        token = strtok(line, " \t\n");

        while ( token != NULL )
        {
            c_token ++;
            p = atof(token);                            // convert string into double
            l->B.dons[c_line][c_token] = p;             // 
            token = strtok(NULL, " \t\n");              // move to next element
        }
    }
    fclose(file);
    if (DEBUG == 1 || DEBUG == 2)     printf("\t\u2713\n");
}

void acceptor_parser(Lambda *l, char *filename) 
{
    if (DEBUG == 1 || DEBUG == 2)     printf("Start getting acceptor site emission Probability:");
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("ERROR: Cannot open file!\n");
        return;
    }
    
    char line[256];
    int row = 0;
    
    while(fgets(line, sizeof(line), file) != NULL)
    {        
        if (line[0] == '%') continue;
                
        char temp[256];
        strcpy(temp, line);
        
        char *token = strtok(temp, " \t\n\r");
        int col = 0;
        
        while(token != NULL && col < 4)
        {
            double val = atof(token);
            l->B.accs[row][col] = val;            
            col++;
            token = strtok(NULL, " \t\n\r");
        }                
        row++;
    }
    
    fclose(file);
    if (DEBUG == 1 || DEBUG == 2)     printf("\t\u2713\n");
}

/* ==================================================== *
 * =============== Emission Probability =============== *
 * ==================================================== */


void exon_intron_parser(Lambda *l, char *filename, int digit)
{
    assert(digit == 0 || digit == 1);                   // 0 for exon, 1 for intron

    if      ( digit == 0 && (DEBUG == 1 || DEBUG == 2) )    printf("Start getting exon   emission  Probability:");
    else if ( digit == 1 && (DEBUG == 1 || DEBUG == 2) )    printf("Start getting intron emission  Probability:");

    FILE *file = fopen(filename, "r");

    char line[256];
    char seq[5];                                        // To hold sequences like AAAA
    double p;

    if (file == NULL)
    {
        if (DEBUG == 1 || DEBUG == 2)     printf("Error: Cannot open file %s\n", filename);
        return;
    }
    
    // Initialize arrays to zero
    if (digit == 0)
    {
        memset(l->B.exon, 0, 256 * sizeof(double));
    } else 
    {
        memset(l->B.intron, 0, 256 * sizeof(double));
    }

    if ( ( fgets(line, sizeof(line), file) != NULL && line[0] == '%' ) && DEBUG == 1)
    {
        printf("\t\u2713");
    }
    else 
    {
        if (DEBUG == 1 || DEBUG == 2)     printf("Warning: No header found in %s\n", filename);
        rewind(file); 
    }

    // Process each line
    while ( fgets( line , sizeof(line), file) != NULL )
    {
        // Skip empty lines or header lines
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '%')   continue;
        
        // Parse the line with sequence and probability
        if (sscanf(line, "%4s %lf", seq, &p) == 2) 
        {
            // Convert sequence to index
            int index = 0;
            for (int i = 0; i < 4; i++) 
            {
                if      (seq[i] == 'A') index = index * 4 + 0;
                else if (seq[i] == 'C') index = index * 4 + 1;
                else if (seq[i] == 'G') index = index * 4 + 2;
                else if (seq[i] == 'T') index = index * 4 + 3;
            }
            
            // Store probability in the appropriate array
            if (index < 256) 
            {
                if      (digit == 0)     l->B.exon[index]   = p;
                else                     l->B.intron[index] = p;
            }
        }
    }
    
    fclose(file);
    if (DEBUG == 1 || DEBUG == 2)     printf("\t\u2713\n");
}

/* ==================================================== *
 * ============== Explicit Duration Prob ============== *
 * ==================================================== */

void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit)
{
    if (digit == 0)
    {
        ed->min_len_exon = -1;
        ed->max_len_exon = 0;
    } else
    {
        ed->min_len_intron = -1;
        ed->max_len_intron = 0;
    }

    FILE *file = fopen(filename, "r");  
    if (!file) { 
        printf("Error: Cannot open file %s\n", filename);
        perror("File open error");
        return; 
    }

    char line[256];
    int duration = 0;
    
    if (DEBUG == 1 || DEBUG == 2) printf("Start parsing explicit duration file: %s (digit=%d)\n", filename, digit);
    
    while (fgets(line, sizeof(line), file))
    {
        // Skip header
        if(line[0] == '%')
        {
            if (DEBUG == 1 || DEBUG == 2)   printf("Skipping header: %s", line);
            continue;
        }
        
        // Parse the probability value
        char *tok = strtok(line, " \t\r\n");
        if(!tok) continue;

        double prob = atof(tok);
        double *arr = (digit == 0) ? ed->exon : ed->intron;

        // Store the probability
        arr[duration] = prob;
        
        // Set minimum length when we encounter first non-zero probability
        if (prob > 0.0)
        {
            if (digit == 0 && ed->min_len_exon == -1)
            {
                ed->min_len_exon = duration;
                if (DEBUG == 1 || DEBUG == 2)   printf("Setting min_len_exon = %d\n", duration);
            }
            if (digit == 1 && ed->min_len_intron == -1)
            {
                ed->min_len_intron = duration;
                if (DEBUG == 1 || DEBUG == 2)   printf("Setting min_len_intron = %d\n", duration);
            }
        }

        duration++;
    }

    // Set maximum length
    if (digit == 0)
    {
        ed->max_len_exon = duration;
        if (DEBUG == 1 || DEBUG == 2)   printf("Set max_len_exon = %d\n", duration);
    } else
    {
        ed->max_len_intron = duration;
        if (DEBUG == 1 || DEBUG == 2)   printf("Set max_len_intron = %d\n", duration);
    }
    
    fclose(file);
    
    if (DEBUG == 1 || DEBUG == 2)   printf("Finished parsing duration file. ");
}