#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

// information we need for isoform creation

typedef struct
{

    int *dons;
    int *accs;
    int *ds_start;
    int *ac_start;
    int *ds_end;
    int *ac_end;
    int dons_count;
    int accs_count;
    int seq_len;
    int flank_size;
    int min_ex;
    int min_in;
    int max_isoform_length;
    int *isoform;

} SpliceSiteData;

/*******************\
prototype being used
\*******************/

// used for analysis seq
void assigner(SpliceSiteData *ssd, const int minin, const int minex, const int flank_size);
void allocate_da_array(SpliceSiteData *ssd, const char *seq);
void splice_site_reader(SpliceSiteData *ssd, const char *seq);

// used for print out array
void pointer_printer(const int *start, const int *end);

// creating isoform
void all_isoform(const SpliceSiteData *ssd, const int *donor, const int *acceptor, const int spot);

// ====================  test sequence ==================== 

char seq[] = "ACGTTGACGTAAGTAAAGCAGCGCCACGAGTAAGAGTAACCGTTTACC";

// ====================  execution     ====================

int main(void){

    SpliceSiteData ssd;

    //remember to change those constant
    //minin, minex, flank
    assigner(&ssd, 2, 2, 2);

    // allocate size to donor and acceptor
    allocate_da_array(&ssd, seq);

    // gathering data from the seq
    splice_site_reader(&ssd, seq);

    // if we wanna print donor and acceptor out for visual
    pointer_printer(ssd.ds_start, ssd.ds_end);
    pointer_printer(ssd.ac_start, ssd.ac_end);

    all_isoform(&ssd, ssd.ds_start, ssd.ac_start, 0);

    free(ssd.dons);
    free(ssd.accs);
    free(ssd.isoform);

    return 0;

}

// set up initial condition about minin, minex, and flank size

void assigner(SpliceSiteData *ssd, const int minin, const int minex, const int flank_size)
{

    ssd->min_in = minin;
    ssd->min_ex = minex;
    ssd->flank_size = flank_size;

}

// asign space to donor and acceptor site based on seq length

void allocate_da_array(SpliceSiteData *ssd, const char *seq)

{
    size_t len = strlen(seq);
    ssd->seq_len = len;
    
    if (len <= 500)

    {
        ssd->dons = malloc ( 20 * sizeof(int) );
        ssd->accs = malloc ( 20 * sizeof(int) );
    }

    else if (len > 500 && len <= 1500)

    {
        ssd->dons = malloc ( 100 * sizeof(int) );
        ssd->accs = malloc ( 100 * sizeof(int) );
    }

    else if (len > 1500 && len <= 3000)

    {
        ssd->dons = malloc ( 150 * sizeof(int) );
        ssd->accs = malloc ( 150 * sizeof(int) );
    }

}

// gt; ag acceptor donor reader

void splice_site_reader(SpliceSiteData *ssd, const char *seq)
{
    ssd->max_isoform_length = floor ( ( ssd->seq_len + 1 - 2 * ssd->flank_size - ssd-> min_ex ) / (ssd->min_ex + ssd->min_in) );
    ssd->isoform = malloc ( ssd->max_isoform_length * sizeof(int) );
    ssd->dons_count = 0;
    ssd->accs_count = 0;


    for (int i = ssd->flank_size+ ssd->min_ex ; i < ssd->seq_len - 1 - ssd->flank_size - ssd->min_ex; i++)
    {

        if (seq[i] == 'G' && seq[i+1] == 'T')
        {
            ssd->dons[ ssd->dons_count ] = i;
            ssd->dons_count++;

        } 

        else if (seq[i] =='A' && seq[i+1] == 'G') 
        {
            ssd->accs[ ssd->accs_count ] = i;
            ssd->accs_count++;
        }
    }

    ssd->ds_start = ssd->dons;
    ssd->ac_start = ssd->accs;
    ssd->ds_end   = ssd->ds_start + ssd->dons_count - 1;
    ssd->ac_end   = ssd->ac_start + ssd->accs_count - 1;

}

/*******************************\
printer session for easier debug
\*******************************/

// print donor and acceptor array with pointer \\

void pointer_printer(const int *start, const int *end)
{

    for (const int *p = start; p <= end; p++ )
    {
        printf("%d\t", *p);
    }

    printf("\n");
    
}

// combinator


void all_isoform(const SpliceSiteData *ssd, const int *donor, const int *acceptor, int spot)
{

    assert(spot % 2 == 0);

    // whenever we are out of donor and acceptor, exit
    if ( donor > ssd->ds_end || acceptor > ssd->ac_end)
    {
        return;
    }

    // creating a nested loop for imitating the intron formation
    for (const int *p1 = donor; p1 <= ssd->ds_end ; p1++ )
    {
        // check the middle exon
        if (spot != 0 && ( ( *p1 - ssd->isoform[spot - 1] ) < (ssd->min_ex + 1) ) )
        {
            continue;
        }

        // ==================================================================\\
        // add whatever features we want to continue the loop if bad isoform \\
        // ==================================================================\\

        // update donor site
        ssd->isoform[spot] = *p1;
        

        // picking acceptors
        for (const int *p2 = acceptor; p2 <= ssd->ac_end; p2++)
        {
            // check the middle intron length
            if ( ( *p2 - *p1 ) < ( ssd->min_in - 1 ) )
            {
                continue;
            }

            // update acceptor site
            ssd->isoform[spot + 1] = *p2;

            pointer_printer( ssd->isoform, ssd->isoform + spot + 1 );

            // ======================================================================\\
            // here we could add the mRNA function which used to store and check then\\
            // ======================================================================\\
            
            // recursion, continue picking
            all_isoform(ssd, p1, p2, spot + 2);
        }
    }   
}


