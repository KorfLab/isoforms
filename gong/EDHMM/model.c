#include "model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

void numerical_transcription(Observed_events *info, const char *seq)
{
    printf("Start transforming original sequence into base4:\n");

    // turns original sequence into int 
    size_t len = strlen(seq);

    info->T = len;
    info->numerical_sequence = malloc ( len * sizeof(int) );

    for( size_t i = 0; i < len; i++)
    {

        // A == 0 , C == 1, G == 2, T == 3  
        if      (seq[i] == 'A')     info->numerical_sequence[i] = 0;
        else if (seq[i] == 'C')     info->numerical_sequence[i] = 1;
        else if (seq[i] == 'G')     info->numerical_sequence[i] = 2;
        else if (seq[i] == 'T')     info->numerical_sequence[i] = 3;
    }

    printf("\tWe get numerical sequence with Seq len: %d\n", info->T);
    printf("\tFinished\n");
    printf("\n");
}

void setup_initial_probability(Lambda *l)                               // actually no longer needed
{
    printf("Start getting initial probability down:");
    l->pi = calloc(HS, sizeof(double) );                                // left-right HMM; only exon are 1
    l->pi[0] = 1;                                                       // initial probability of exon are 1
    printf("\t\u2713\n");
}

int power(int base, int exp)                                            // wtf, C don't have power for int
{
    int result = 1;

    for ( int i = 0 ; i < exp ; i++ )
    {
        result *= base;
    }
    return result;
}

int base4_to_int(int *array, int beg, int length) 
{
    int values = 0;
    int value;

    for (int i = 0; i < length; i++)
    {
        value  =  array[beg + i];
        values += value * power(4, length - i - 1);
    }
    
    return values;
}

double total_prob(double *array, int length)
{
    double value = 0.0;

    for (int i = 0 ; i < length ; i ++)
    {
        if (array[i] == 0.0)
        {
            value = 0.0;
            break;
        }

        value += log(array[i]);
    }

    if (value != 0.0)   value = exp(value);

    return value;
}

void initialize_donor_transition_matrix(Lambda *l, Apc *a, int depth)   // set the depth to 0 initially
{   
    if (depth == 5)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 5);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 5);                          // get total prob
        l->A.dons[index]  = value;                                      // store the value
        return;
    }

    for ( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.dons[depth][i];                                 // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_donor_transition_matrix(l, a, depth + 1);            // send into next node
    }
}

void initialize_acceptor_transition_matrix(Lambda *l, Apc *a, int depth)// set the depth to 0 initially
{   
    if (depth == 6)                                                     // exit recursion; calculation start
    {
        int index = base4_to_int(a->position, 0, 6);                    // this is where we plan to store that value
        double value = total_prob(a->prob, 6);                          // get total prob
        l->A.accs[index]  = value;                                      // store the value
        return;
    }

    for ( int i = 0; i < 4 ; i++ )
    {
        double p = l->B.accs[depth][i];                                 // get the exon emission prob for current node
        a->prob[depth] = p;                                             // record the emission prob inside the array
        a->position[depth] = i;                                         // record each base pair each node choose
        initialize_acceptor_transition_matrix(l, a, depth + 1);         // send into next node
    }
}

void normalize_transition_prob(Lambda *l, int len, int dons_or_accs)
{
    printf("Start normalizing transition prob:");

    double sum = 0.0;

    if (dons_or_accs == 0)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.dons[i];

        double value;
        for ( int i = 0 ; i < len ; i ++ )
        {
            value = l->A.dons[i];
            if (value == 0.0)       continue;
            else                    l->A.dons[i] = exp ( log(value) - log(sum) );
            
        }
    }
    else if (dons_or_accs == 1)
    {
        for ( int i = 0 ; i < len ; i ++ )  sum += l->A.accs[i];

        double value;
        for ( int i = 0 ; i < len ; i ++ )
        {
            value = l->A.accs[i];

            if (value == 0.0)       continue;
            else                    l->A.accs[i] = exp( (log(value) - log(sum) ) );
        }
    }

    printf("\tFinished\n");
}

double log_sum_exp(double *array, int n) 
{
    double max = array[0];

    for (int i = 1; i < n ; i++)
    {
        if( array[i] > max )     max = array[i];                  
    }

    double sum = 0.0;

    for( int i = 0 ; i < n ; i++)
    {
        sum += exp(array[i] - max);
    }

    return max + log(sum);       
}

void allocate_alpha(Observed_events *info, Forward_algorithm *alpha , Explicit_duration *ed)                            
{
    printf("Start allocate memory for the forward algorithm:");

    int arary_size = info->T - 2 * FLANK;

    /*
        alpha->a[t][i]
        [t]: specific time among all observed events 
        [i]: [0] for exon ; [1] for intron

        each spot is storing a(t)(m, 1) ; based on 2006 implementation
        [m]: types of hidden state
    */

    alpha->a    = malloc ( ( arary_size ) * sizeof(double*) ); 

    for (int i = 0 ; i < arary_size; i++ )     
        alpha->a[i] = calloc( HS , sizeof(double) );                                        
    
    /*
        alpha->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
        based on forward algorithm; each at(m, d) partially depends one a(t-1)(m, d+1)
    */

    alpha->basis    = malloc( HS * sizeof(double*) );                   
    alpha->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    alpha->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    printf("\tFinished\n");
}

void basis_forward_algorithm(Lambda *l, Explicit_duration *ed,  Forward_algorithm *alpha, Observed_events *info)
{
    printf("Start forward algorithm basis calculation:");

    /*
        first loop: updating all possible a(1)(m, d)

        a(1)(m, d) = pi(m) * bm(o1) * pm(d)
        pi(m):  initial probability
        bm(o1): emission probability
        pm(d): explicit duration probability
        [tau]: used for stimulate different duration

        well we don't need initial prob for us; since uh, exon = 1.0; intron = 0.0;
    */

    int tau_exon;
    int comp_len_exon = info->T - 2 * FLANK;

    if  (comp_len_exon > ed->max_len_exon)      tau_exon = ed->max_len_exon;
    else                                        tau_exon = comp_len_exon;

    int    index         = base4_to_int(info->numerical_sequence, FLANK - 3, 4);
    double emission_prob = l->B.exon[index];

    /*
        explicit duration is the only term that gonna be 0.0
        if it's 0; means it's invalid length for exon
        directly assign it as 0.0
    */

    double total;
    double ed_prob;

    for( int d = 0 ; d < tau_exon ; d ++)
    {
        ed_prob = ed->exon[d];

        if (ed_prob == 0.0)
        {
            alpha->basis[0][d] = 0.0;
        }
        else
        {
            total = exp( log(ed_prob) + log(emission_prob) );
            alpha->basis[0][d] = total;
        }
    }

    /*
        this part is assign prob for intron
        since for this model; it;s kinda different
        we don't have initial intron probability
        btw the real intron initial probability comes after min len of exon
        which make sense
        so we make our first layer of calculation at point where passes min exon len
    */

    int tau_intron;
    int comp_len_intron = comp_len_exon - ed->min_len_exon * 2;

    if ( comp_len_intron > ed->max_len_intron)      tau_intron = ed->max_len_intron;
    else                                            tau_intron = comp_len_intron;

    index = base4_to_int(info->numerical_sequence, FLANK + ed->min_len_exon - 3 , 4);
    emission_prob = l->B.intron[index];

    for( int d = 0 ; d < tau_intron ; d ++)
    {
        ed_prob = ed->intron[d];

        if (ed_prob == 0.0)
        {
            alpha->basis[1][d] = 0.0;
        }
        else
        {
            total = exp( log(ed_prob) + log(emission_prob) );
            alpha->basis[1][d] = total;
        }
    }

    printf("\tFinished\n");
}

void forward_algorithm(Lambda *l, Forward_algorithm *alpha, Observed_events *info, Explicit_duration *ed)
{
    printf("Start computation for forward algorithm:");

    /*
        recall
        [tau]: the residual/remaining time for explicit duration
    */

    int len = info->T - 2 * FLANK;
    int tau = len;

    for ( int t = 1 ; t < info->T - 2 * FLANK ; t ++ )
    {
        tau --;

        for ( int i = 0 ; i < HS ; i ++ )
        {
            int tau_modified;
            if ( i == 1)    tau_modified = tau - ed->min_len_exon;
            else            tau_modified = tau;

            /*
                boundary condition

                [t < ed->min_len_exon]: before first min_len_exon; no intron avaliable
                [t > len - ed->min_len_exon]: no intron are possible unless a min_len_exon able to exist in the end

                for computation
                
                we just need to record those α(t)(intron, 1)
                we don't need to update the value in the alpha->basis ; the layer of network
            */

            if ( ( t < ed->min_len_exon || t > len - ed->min_len_exon ) && i == 1)
            {
                alpha->a[t][i] = 0.0;
                continue;
            }

            double trans_prob;
            int index_trans_prob;

            double node_trans;
            double alpha_trans;
            double ed_prob;
            double node_continue;

            double emission_prob;
            int index_emission_prob;

            int j = (i == 0) ? 1 : 0;
            double total;

            /*
                first part
                α(t)(m, d) = bm(ot) * ( α(t - 1)(m, d + 1) + α(t - 1)(n, 1) * a(nm) * pm(d))
                    for d > 1

                [trans_prob]: a(nm)
                [node_trans]: α(t - 1)(n , 1)
                [alpha_trans]: α(t - 1)(n , 1) * a(nm)      aka: trans_prob * node_trans
                [ed_prob]:  pm(d)
                [node_continue]: α(t - 1)(m, d + 1)
                [emission_prob]: bm(ot)
                [j]: conjudated hidden state                aka: i = exon; j = intron | i = intron; j = exon
                [total]: everything without bm(ot)          aka: α(t - 1)(m, d + 1) + α(t - 1)(n, 1) * a(nm) * pm(d)
            */

            if ( i == 0 )
            {
               index_trans_prob = base4_to_int(info->numerical_sequence , t + FLANK - 6, 6);
               trans_prob = l->A.accs[index_trans_prob];
            }
            else
            {
               index_trans_prob = base4_to_int(info->numerical_sequence , t + FLANK , 5);
               trans_prob = l->A.dons[index_trans_prob];
            }

            node_trans = alpha->a[t - 1][j];

            if      (node_trans == 0.0)           alpha_trans = 0.0;
            else if (trans_prob == 0.0)           alpha_trans = 0.0;
            else                                  alpha_trans = exp( log(trans_prob) + log(node_trans) );
            
            index_emission_prob = base4_to_int(info->numerical_sequence, t + FLANK - 3, 4);
            emission_prob = ( i == 0 ) ? l->B.exon[index_emission_prob] : l->B.intron[index_emission_prob];

            for ( int d = 1 ; d < tau_modified ; d ++ )
            {   
                node_continue = alpha->basis[i][d + 1];

                if ( node_continue == 0.0 )     l->log_values[0] = 0.0;
                else                            l->log_values[0] = node_continue;

                ed_prob = ( i == 0 ) ? ed->exon[d] : ed->intron[d];

                if   (ed_prob == 0.0)           l->log_values[1] = 0.0;
                else                            l->log_values[1] = exp ( log(alpha_trans) + log(ed_prob) );

                total = log_sum_exp(l->log_values, 2);
                alpha->basis[i][d] = exp( log(total) + log(emission_prob) );
            }

            /*
                second part
                α(t)(m, 1) = α(t - 1)(m, 2) * bm(ot)
                    since mostly explicit duration when d = 0 is 0 ; the rest terms canceled
                    and α(t)(m, 1) for both intron and exon is stored in alpha->basis[t][hidden state]

                [node_continue]: α(t - 1)(m, 2)     reused variable here
            */

            node_continue = alpha->basis[i][1];

            if ( node_continue  == 0.0 )    total = 0.0;
            else                            total = exp( log( node_continue ) + log(emission_prob) );

            alpha->a[t][i]     = total;
            alpha->basis[i][0] = total;
        }
    }
    printf("\tFinished\n");
}

void free_alpha(Observed_events *info, Forward_algorithm *alpha)
{
    printf("Clearing up forward algorithm memory:");
    
    int array_size = info->T - 2 * FLANK;
    
    for (int i = 0; i < array_size; i++)        free(alpha->a[i]);
    
    free(alpha->a);
    free(alpha->basis[0]);
    free(alpha->basis[1]);
    free(alpha->basis);

    printf("\tFinished\n");
}

void allocate_viterbi(Viterbi_algorithm *vit, Observed_events *info)
{
    printf("Start Initialize Viterbi Algorithm");

    /*
        recursive viterbi formula
        γ(t) = y(t+1)(m) + sum(n != m) ( ξ(t+1)(m, n) - ξ(t+1)(n, n) )

        for our model: degrade that formula
        γ(t) = y(t+1)(m) + ξ(t+1)(m, n) - ξ(t+1)(n, n)

        given definition of ξ
        ξ(t)(m, n) = α(t - 1)(m, 1) * a(mn) * bn(Ot) * sum(d >= 1) (pn(d) β(t)(n, d))

        in terms of t + 1
        ξ(t+1)(m, n) = α(t)(m, 1) * a(mn) * bn(Ot + 1) * sum(d >= 1) (pn(d) β(t+1)(n, d))
    */

    vit->xi    = calloc( HS , sizeof(double) );
    vit->gamma = malloc( HS * sizeof(double) );
    vit->path  = malloc( (info->T - 2 * FLANK) * sizeof(int) );

    int array_size = info->T - 2 * FLANK;

    vit->xi_sum= malloc ( HS * sizeof(double*) ); 

    for (int i = 0 ; i < HS; i++ )     
        vit->xi_sum[i] = calloc( array_size , sizeof(double) );      

    printf("\tFinished\n");
}

void argmax_viterbi(Viterbi_algorithm *vit, int t)
{   
    int argmax;

    /*
        for our model: consider following viterbi formula

        γ(t)(exon)   = γ(t+1)(exon)   + ξ(exon, intron) - ξ(intron, exon)
        γ(t)(intron) = γ(t+1)(intron) + ξ(intron, exon) - ξ(exon, intron)
        
        deduct from 2006 paper
    */

    vit->gamma[0] += vit->xi[0] - vit->xi[1];
    vit->gamma[1] += vit->xi[1] - vit->xi[0];

    vit->xi_sum[0][t] = vit->xi[0];
    vit->xi_sum[1][t] = vit->xi[1];

    if      ( vit->gamma[0] > vit->gamma[1] )    argmax = 0;
    else if ( vit->gamma[0] < vit->gamma[1] )    argmax = 1;
    else printf("\nDoes this really gonna happen? At %d. γ[0]: %f γ[1]: %f", t , vit->gamma[0], vit->gamma[1]);

    vit->path[t] = argmax;
}

void xi_calculation(Lambda *l, Forward_algorithm *alpha, Viterbi_algorithm *vit, Observed_events *info, double backward_sum, int t, int type)
{
    assert(type == 0 || type == 1);

    /*
        input parameter

        [backward_sum]: (sum d>= 1)[ pn(d) * β(n, d)]
            compute inside each layer of backward algorithm

        [type]: exon 0 or intron 1
            so we know which ξ(m, n) we shall update in vit
        
    */

    double alpha_component;
    
    double trans_prob;
    int    index_trans_prob;

    double emission_prob;
    int    index_emission_prob;

    double xi;

    /*
        formula
        ξ(t)(m, n) = α(t - 1)(m, 1) * a(mn) * bn(ot) * (sum d>= 1)[ pn(d) * β(n, d)]   

        update formula
        ξ(t+1)(m, n) = α(t)(m, 1) * a(mn) * bn(o t+1) * (sum d>= 1)[ pn(d) * β(n, d)] 

        [alpha_component]: α(t)(m, 1)
        [trans_prob]: a(mn)
        [emission_prob]: bn(o t+1)
    */

    alpha_component = alpha->a[t - 1][type];

    if (type == 0)
    {
        index_trans_prob = base4_to_int(info->numerical_sequence, t + FLANK, 5);
        trans_prob = l->A.dons[index_trans_prob];
    }
    else
    {
        index_trans_prob = base4_to_int(info->numerical_sequence, t + FLANK - 6, 6);
        trans_prob = l->A.accs[index_trans_prob];
    }
        
    index_emission_prob = base4_to_int(info->numerical_sequence, t + FLANK - 3, 4);
    emission_prob = (type == 0) ? l->B.intron[index_emission_prob] : l->B.exon[index_emission_prob];

    if      (trans_prob == 0.0)         xi = 0.0;
    else if (alpha_component == 0.0)    xi = 0.0;
    else if (backward_sum == 0.0)       xi = 0.0;
    else    xi = exp( log(trans_prob) + log(alpha_component) + log(emission_prob) + log(backward_sum) );

    vit->xi[type] = xi;
}

void viterbi_basis(Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    printf("Start assign basis for Viterbi Algorithm:");

    /*
        γ(t)(m) = sum(d>=1) α(t)(m, d)
        at final t; there is only one α(t)(m, 1) existed
        based on definition for d
        (if forward algorithm is compute without boundary issue)

        γ(t)(1) = shall be 0; intron uh is invalid in that position
    */

    double gamma_exon;
    double gamma_intron;

    gamma_exon   = alpha->basis[0][0];
    gamma_intron = 0.0;

    vit->gamma[0] = gamma_exon;
    vit->gamma[1] = gamma_intron;
    
    printf("\tFinished\n");
}

void allocate_beta(Backward_algorithm *beta, Explicit_duration *ed)                             
{
    printf("Start allocate memory for the backward algorithm:");
                                    
    /*
        β->basis[i][d]
        [i]: [0] for exon ; [1] for intron
        [d]: max duration for exon or intron ; depends on [i]

        each one assign 1D array for each t - 1 layer of all possible D computation
    */

    beta->basis    = malloc( HS * sizeof(double*) );                   
    beta->basis[0] = calloc( ed->max_len_exon, sizeof(double) );
    beta->basis[1] = calloc( ed->max_len_intron, sizeof(double));

    printf("\tFinished\n");
}

void initial_backward_algorithm(Backward_algorithm *beta)
{
    printf("Start initialize backward algorithm:");

    for ( int i = 0 ; i < HS ; i++ )
    {
        beta->basis[i][0] = 1.0;
    }

    printf("\tFinished\n");
}

void backward_algorithm(Lambda *l, Backward_algorithm *beta, Observed_events *info, Explicit_duration *ed, Viterbi_algorithm *vit, Forward_algorithm *alpha)
{
    printf("Start Backward Algorithm:");

    int start_bps = info->T - 2 * FLANK - 1;
    int tau = 0;

    printf("\n\tStart at first base pair location: %d\n", start_bps);

    for ( int t = start_bps ; t >= 0 ; t-- )                                      // -1 cuz array start at 0; -1 again since already set up last one
    {
        argmax_viterbi(vit, t);
        
        if (t == 0)     break;

        tau ++;

        for ( int i = 0 ; i < HS ; i ++ )                                                       // the before position; 0 for exon, 1 for intron
        {
            int tau_modified;

            if ( i == 0 )   tau_modified = tau;
            else            tau_modified = tau - ed->min_len_exon;

            if       ( i == 0 && tau >= ed->max_len_exon)                        tau_modified = ed->max_len_exon;
            else if  ( i == 1 && tau >= ed->max_len_intron - ed->min_len_exon)   tau_modified = ed->max_len_intron;
    
            /*
                boundary condition

                [final bps]: total_bps - 2 * FLANK (- 1) as in array form
                [last intron range]: final bps - min_exon (-1) in array form

                [first intron range]: FLANK + min_exon

                in terms of backward algorithm; purpose is calculate summation of all possible backward transition
                with in those boundary; directly assign such calculation to 0; represent they are inreasonable nodes for network
            */

            double backward_sum;

            if ( ( tau_modified <= 0 || t < FLANK + ed->min_len_exon ) && i == 1) 
            {
                backward_sum = 0.0;
                xi_calculation(l, alpha, vit, info, backward_sum, t, 1);
                continue;
            }

            double ed_prob;
            double possible_node;

            double trans_prob;
            int    index_trans_prob;

            double emission_prob;
            int    index_emission_prob;

            double total;

            /*
                formula
                β(t)(m, 1) = sum(n != m) * a(mn) * bn(Ot+1) * sum(d>=1) pn(d) * β(t+1)(n, d)

                [backward_sum]: sum(d>=1) pn(d) * β(t + 1)(n , d)
                    can be coupled with ξ calculation for viterbi algorithm

                [ed_prob]: pn(d)
                [possible_node]: β(t + 1)(n , d)
                [trans_prob]: a(mn)
                [emission_prob]: bn(ot+1)
                [j]: conjugated hidden state
                    if exon(0) it's intron(1); vice versa
            */

            int j = ( i == 0) ? 1 : 0;

            for ( int d = 0 ; d < tau_modified ; d++ )
            {
                ed_prob = ( j == 0 ) ? ed->exon[d] : ed->intron[d];
                possible_node = beta->basis[j][d];

                if (ed_prob == 0.0)     l->log_values[d] = 0.0;
                else                    l->log_values[d] = exp( log( possible_node ) + log (ed_prob) );
            }

            backward_sum = log_sum_exp(l->log_values, tau_modified);
            xi_calculation(l, alpha, vit, info, backward_sum, t, j);

            if ( i == 0 )
            {   
                index_trans_prob = base4_to_int(info->numerical_sequence, t + 1 + FLANK , 5);
                trans_prob = l->A.dons[index_trans_prob];
            }
            else
            {
                index_trans_prob = base4_to_int(info->numerical_sequence, t - 5 + FLANK , 6);
                trans_prob = l->A.accs[index_trans_prob];
            }

            index_emission_prob = base4_to_int(info->numerical_sequence, t - 2 + FLANK, 4);
            emission_prob = (j == 0) ? l->B.exon[index_emission_prob] : l->B.intron[index_emission_prob];

            if   (trans_prob == 0.0)   total = 0.0;
            if   (backward_sum == 0.0) total = 0.0;
            else                       total = exp( log(trans_prob) + log(emission_prob) + log(backward_sum) );

            /*
                β(t)(m, d) = bm(Ot+1) * β(t+1)(m, d - 1)
                    for all possible d > 1

                [emission_prob]: bm(Ot+1)       notice: this is not conjugated state right here
            */
            
            index_emission_prob = base4_to_int(info->numerical_sequence, t - 2 + FLANK, 4);
            emission_prob = (i == 0) ? l->B.exon[index_emission_prob] : l->B.intron[index_emission_prob];

            /*
                logic here is linear array for all layer of network
                and each time passed there would be one node lost for each layer
                every computation based on the previous node
                so we cannot directly update that after calculation

                [first_node]: β(t)(m, 1)
                    it is always the first node calculated
                [previous_node]: record β(t)(m, d - 1) for next β(t)(m, d) calculation before update them
            */

            double first_node;
            double previous_node;

            first_node = total;
            previous_node = total;

            for( int d = 1 ; d < tau_modified ; d++ )
            {   

                if      (previous_node == 0.0)  total = 0.0;
                else                            total = exp( log(emission_prob) + log(previous_node) );

                previous_node = beta->basis[i][d];
                beta->basis[i][d] = total;
            } 

            beta->basis[i][0] = first_node;
        }
    }

    vit->xi_sum_exon   = log_sum_exp(vit->xi_sum[0], info->T - 2 * FLANK);
    vit->xi_sum_intron = log_sum_exp(vit->xi_sum[1], info->T - 2 * FLANK);

    printf("\tThis is xi sum for exon throughout the time %f\n",   vit->xi_sum_exon);
    printf("\tThis is xi sum for intron throughout the time %f\n", vit->xi_sum_intron);
    printf("\tFinished.\n");
}

void free_beta(Backward_algorithm *beta)
{
    printf("Clearning up backward algorithm memory:");

    free(beta->basis[0]);
    free(beta->basis[1]);
    free(beta->basis);

    printf("\tFinished\n");
}

void free_viterbi(Viterbi_algorithm *vit)
{
    printf("Clearning up viterbi algorithm memory:");

    free(vit->path);
    free(vit->gamma);
    free(vit->xi);
    free(vit->xi_sum[0]);
    free(vit->xi_sum[1]);
    free(vit->xi_sum);

    printf("\tFinished\n");
}

void viterbi_path_test(Viterbi_algorithm *vit, Observed_events *info, Explicit_duration *ed)
{
    printf("\nStart Viterbi Check:\n");
    
    int segment = 1;
    int start_pos = FLANK;
    int current_state = vit->path[0];
    int transition_pos;
    int length;
    
    int exon_count = 0;
    int intron_count = 0;
    int invalid_exons = 0;
    int invalid_introns = 0;
    int canonical_donors = 0;
    int canonical_acceptors = 0;
    int total_donors = 0;
    int total_acceptors = 0;
    
    printf("%-8s %-8s %-8s %-8s %-8s %-12s %-25s\n", 
           "Segment", "Type", "Start", "End", "Length", "Valid", "Splice Site");
    printf("%-8s %-8s %-8s %-8s %-8s %-12s %-25s\n", 
           "-------", "----", "-----", "---", "------", "--------", "-----------");
    
    if (current_state == 0) 
    {
        exon_count++;
        printf("%-8d %-8s %-8d ", segment, "Exon", start_pos);
    } else 
    {
        intron_count++;
        printf("%-8d %-8s %-8d ", segment, "Intron", start_pos);
    }
    
    for (int i = 1; i < info->T - 2 * FLANK; i++)
    {
        if (vit->path[i] != current_state)
        {
            transition_pos = i + FLANK;
            length = transition_pos - start_pos;
            
            int valid_len = 0;

            if (current_state == 0)
            { 
                valid_len = (length >= ed->min_len_exon);
                if (!valid_len) invalid_exons++;
            } else 
            {
                valid_len = (length >= ed->min_len_intron);
                if (!valid_len) invalid_introns++;
            }
            
            printf("%-8d %-8d %-12s ", transition_pos - 1, length, valid_len ? "Yes" : "No");
            
            if (current_state == 0) 
            {
                // Exon to Intron transition (donor site)
                total_donors++;
                
                // The first base of the intron should be at transition_pos
                // We want to show transition_pos and the next 4 bases (total 5 bases)
                if (transition_pos + 4 < info->T)
                {
                    char donor_site[6];
                    // FIXED: Start exactly at the intron start (no off-by-one)
                    strncpy(donor_site, &info->original_sequence[transition_pos], 5);
                    donor_site[5] = '\0';
                    printf("Donor: %s ", donor_site);
                    
                    // Check for canonical GT at the first two positions of the intron
                    if (info->original_sequence[transition_pos] == 'G' && 
                        info->original_sequence[transition_pos + 1] == 'T') {
                        printf("(Canonical GT)\n");
                        canonical_donors++;
                    } else {
                        printf("(Non-canonical)\n");
                    }
                } else {
                    printf("Donor: (insufficient seq)\n");
                }
            } else { 
                // Intron to Exon transition (acceptor site)
                total_acceptors++;
                
                // We need the last 6 bases of the intron (ending at transition_pos-1)
                if (transition_pos >= 6) {
                    char acceptor_site[7];
                    strncpy(acceptor_site, &info->original_sequence[transition_pos - 6], 6);
                    acceptor_site[6] = '\0';
                    printf("Acceptor: %s ", acceptor_site);
                    
                    // Check for canonical AG at the last two positions of the intron
                    if (info->original_sequence[transition_pos - 2] == 'A' && 
                        info->original_sequence[transition_pos - 1] == 'G') {
                        printf("(Canonical AG)\n");
                        canonical_acceptors++;
                    } else {
                        printf("(Non-canonical)\n");
                    }
                } else {
                    printf("Acceptor: (insufficient seq)\n");
                }
            }
            
            // Update for next segment
            current_state = vit->path[i];
            start_pos = transition_pos;
            segment++;
            
            // Count new segment
            if (current_state == 0) {
                exon_count++;
                printf("%-8d %-8s %-8d ", segment, "Exon", start_pos);
            } else {
                intron_count++;
                printf("%-8d %-8s %-8d ", segment, "Intron", start_pos);
            }
        }
    }
    
    // Process final segment
    int end_pos = info->T - FLANK;
    length = end_pos - start_pos;
    
    // Check segment validity
    int valid_len = 0;
    if (current_state == 0) { // Exon
        valid_len = (length >= ed->min_len_exon);
        if (!valid_len) invalid_exons++;
    } else { // Intron
        valid_len = (length >= ed->min_len_intron);
        if (!valid_len) invalid_introns++;
    }
    
    printf("%-8d %-8d %-12s (end of sequence)\n", 
          end_pos - 1, length, valid_len ? "Yes" : "No");
    
    // Print summary
    printf("\nViterbi Path Summary:\n");
    printf("Total segments: %d\n", segment);
    printf("Exons: %d (Invalid length: %d)\n", exon_count, invalid_exons);
    printf("Introns: %d (Invalid length: %d)\n", intron_count, invalid_introns);
    printf("Canonical donor sites (GT): %d/%d\n", canonical_donors, total_donors);
    printf("Canonical acceptor sites (AG): %d/%d\n", canonical_acceptors, total_acceptors);
    
    // Check gene structure
    if (vit->path[0] != 0) {
        printf("\nWARNING: Prediction starts with an intron (biologically unusual)!\n");
    }
    if (current_state != 0) {
        printf("\nWARNING: Prediction ends with an intron (biologically unusual)!\n");
    }
}