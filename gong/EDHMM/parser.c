#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "model.h"

// read the sequence from file

void read_sequence_file(const char *filename, Observed_events *info)
{
    printf("Start reading the sequence data:\n");

    FILE *file = fopen(filename, "r");
    
    if (file == NULL)
    {
        printf("Error: Cannot open sequence file %s\n", filename);
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
    
    // Parse each line
    char *line = strtok(buffer, "\n");
    while (line != NULL)
    {
        
        // Extract valid DNA characters
        for (int i = 0; line[i] != '\0'; i++)
        {
            if (line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T')
            {
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

    printf("\tWe get original sequence with Seq len: %zu\n", seq_index);
    printf("\tFinished\n");
    printf("\n");
}

// emission probability //

void donor_parser(Lambda *l, char *filename)            // get emission probability for donor site
{
    printf("Start getting donor site emission Probability:");
    FILE *file = fopen(filename, "r");

    char line[256];
    char *token;
    double p;                                           // probability we are going to store

    if (file == NULL)
    {
        printf("Can't find file for donor site emission probability!\n");
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
    printf("\t\u2713\n");
}

void acceptor_parser(Lambda *l, char *filename)         // get emission probability for acceptor site
{
    printf("Start getting acceptor site emission Probability:");

    FILE *file = fopen(filename, "r");

    char line[256];
    char *token;
    double p;                                           // probability we are gonna store

    if (file == NULL)
    { 
        printf("Can't find file for donor site emission probability!\n");
        return;
    }
    
    int c_line = -1;                                    // count of line

    while( fgets( line, sizeof(line) , file) != NULL )  // nest while loop to get elements
    {
        if ( line[0] == '%')     continue;              // skip the first line   

        c_line++;
        int c_token = -1;

        token = strtok(line, " \t\n");                  // get each probability

        while ( token != NULL )                         
        {
            c_token ++;
            p = atof(token);                            // convert string into double
            l->B.accs[c_line][c_token] = p;             // store the value
            token = strtok(NULL, " \t\n");              // move to next element
        }
    }
    fclose(file);

    printf("\t\u2713\n");
}

void exon_intron_parser(Lambda *l, char *filename, int digit)
{
    assert(digit == 0 || digit == 1);                   // 0 for exon, 1 for intron

    if      ( digit == 0 )    printf("Start getting exon   emission  Probability:");
    else if ( digit == 1 )    printf("Start getting intron emission  Probability:");

    FILE *file = fopen(filename, "r");

    char line[256];
    char seq[5];                                        // To hold sequences like AAAA
    double p;

    if (file == NULL)
    {
        printf("Error: Cannot open file %s\n", filename);
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

    // Skip header line
    if (fgets(line, sizeof(line), file) != NULL && line[0] == '%')      printf("\t\u2713");
    else 
    {
        printf("Warning: No header found in %s\n", filename);
        rewind(file); // Go back to start if no header
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
    printf("\t\u2713\n");
}

// eplicit_duration //

void explicit_duration_probability(Explicit_duration *ed, char *filename, int digit)
{
    assert(digit == 0 || digit == 1);                   // 0 for exon, 1 for intron

    if      ( digit == 0 ) printf("Starting getting exon explicit duration probability");
    else if ( digit == 1 ) printf("Starting getting intron explicit duration probability");
    
    FILE *file = fopen(filename, "r");


    if (file == NULL)
    {
        printf("Error opening explicit duration file: %s\n", filename);
        perror("Error details");
        return;
    }

    char line[256];
    char *token;
    double p;
    
    // Initialize arrays to zero
    if (digit == 0) 
    {
        memset(ed->exon, 0, 1000 * sizeof(double));
        ed->min_len_exon = 0;
        ed->max_len_exon = 0;
    } else 
    {
        memset(ed->intron, 0, 1000 * sizeof(double));
        ed->min_len_intron = 0;
        ed->max_len_intron = 0;
    }

    int c_line = 0;
    int nonzero_values = 0;

    while (fgets(line, sizeof(line), file) != NULL)
    {
        if (line[0] == '%')     continue;  // Skip header line

        c_line++;
        
        // Clean up the line - remove whitespace
        token = strtok(line, " \t\r\n");
        if (token == NULL) continue;
        
        p = atof(token);
        
        // Store the probability and track first non-zero value
        if (digit == 0) 
        {
            ed->exon[c_line] = p;
            if (p > 0 && ed->min_len_exon == 0)     ed->min_len_exon = c_line;
        } else 
        {
            ed->intron[c_line] = p;
            if (p > 0 && ed->min_len_intron == 0)   ed->min_len_intron = c_line;
        }
        
        if (p > 0) nonzero_values++;
    }

    fclose(file);
    printf("\t\u2713\n");

    // Set max length
    if (digit == 0) 
    {
        ed->max_len_exon = c_line;
        printf("\tExon duration: min=%d, max=%d, found %d non-zero values\n\n", 
               ed->min_len_exon, ed->max_len_exon, nonzero_values);
    } else 
    {
        ed->max_len_intron = c_line;
        printf("\tIntron duration: min=%d, max=%d, found %d non-zero values\n\n", 
               ed->min_len_intron, ed->max_len_intron, nonzero_values);
    }

}