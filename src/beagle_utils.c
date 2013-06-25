/*
 * author: Imran Fanaswala
 */

#ifndef BEAGLE_UTILS_CPP
#define BEAGLE_UTILS_CPP

#include  <stdio.h>
#include "beagle_utils.h"

double* tips_to_partials_nucl(const char *sequence) {
    int n = strlen(sequence);
    double *partials = (double*)malloc(sizeof(double) * n * 4);

    int k = 0;
    for (int i = 0; i < n; i++) {
        switch (sequence[i]) {
            case 'A':
                partials[k++] = 1;
                partials[k++] = 0;
                partials[k++] = 0;
                partials[k++] = 0;
                break;
            case 'C':
                partials[k++] = 0;
                partials[k++] = 1;
                partials[k++] = 0;
                partials[k++] = 0;
                break;
            case 'G':
                partials[k++] = 0;
                partials[k++] = 0;
                partials[k++] = 1;
                partials[k++] = 0;
                break;
            case 'T':
                partials[k++] = 0;
                partials[k++] = 0;
                partials[k++] = 0;
                partials[k++] = 1;
                break;
            case 'X':
                partials[k++] = 1;
                partials[k++] = 1;
                partials[k++] = 1;
                partials[k++] = 1;
                break;
            default:
                Warn_And_Exit("TODO: Add other states");
        }
    }
    return partials;
}

void print_beagle_flags(long inFlags) {
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CPU)      fprintf(stdout, " PROCESSOR_CPU");
    if (inFlags & BEAGLE_FLAG_PROCESSOR_GPU)      fprintf(stdout, " PROCESSOR_GPU");
    if (inFlags & BEAGLE_FLAG_PROCESSOR_FPGA)     fprintf(stdout, " PROCESSOR_FPGA");
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CELL)     fprintf(stdout, " PROCESSOR_CELL");
    if (inFlags & BEAGLE_FLAG_PRECISION_DOUBLE)   fprintf(stdout, " PRECISION_DOUBLE");
    if (inFlags & BEAGLE_FLAG_PRECISION_SINGLE)   fprintf(stdout, " PRECISION_SINGLE");
    if (inFlags & BEAGLE_FLAG_COMPUTATION_ASYNCH) fprintf(stdout, " COMPUTATION_ASYNCH");
    if (inFlags & BEAGLE_FLAG_COMPUTATION_SYNCH)  fprintf(stdout, " COMPUTATION_SYNCH");
    if (inFlags & BEAGLE_FLAG_EIGEN_REAL)         fprintf(stdout, " EIGEN_REAL");
    if (inFlags & BEAGLE_FLAG_EIGEN_COMPLEX)      fprintf(stdout, " EIGEN_COMPLEX");
    if (inFlags & BEAGLE_FLAG_SCALING_MANUAL)     fprintf(stdout, " SCALING_MANUAL");
    if (inFlags & BEAGLE_FLAG_SCALING_AUTO)       fprintf(stdout, " SCALING_AUTO");
    if (inFlags & BEAGLE_FLAG_SCALING_ALWAYS)     fprintf(stdout, " SCALING_ALWAYS");
    if (inFlags & BEAGLE_FLAG_SCALING_DYNAMIC)    fprintf(stdout, " SCALING_DYNAMIC");
    if (inFlags & BEAGLE_FLAG_SCALERS_RAW)        fprintf(stdout, " SCALERS_RAW");
    if (inFlags & BEAGLE_FLAG_SCALERS_LOG)        fprintf(stdout, " SCALERS_LOG");
    if (inFlags & BEAGLE_FLAG_VECTOR_NONE)        fprintf(stdout, " VECTOR_NONE");
    if (inFlags & BEAGLE_FLAG_VECTOR_SSE)         fprintf(stdout, " VECTOR_SSE");
    if (inFlags & BEAGLE_FLAG_VECTOR_AVX)         fprintf(stdout, " VECTOR_AVX");
    if (inFlags & BEAGLE_FLAG_THREADING_NONE)     fprintf(stdout, " THREADING_NONE");
    if (inFlags & BEAGLE_FLAG_THREADING_OPENMP)   fprintf(stdout, " THREADING_OPENMP");
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_CPU)      fprintf(stdout, " FRAMEWORK_CPU");
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_CUDA)     fprintf(stdout, " FRAMEWORK_CUDA");
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_OPENCL)   fprintf(stdout, " FRAMEWORK_OPENCL");
    fflush(stdout);
}

void print_beagle_resource_list()
{
    BeagleResourceList* rList;
    rList = beagleGetResourceList();
    fprintf(stdout, "Available resources:\n");
    for (int i = 0; i < rList->length; i++) {
        fprintf(stdout, "\tResource %i:\n\t\tName : %s\n", i, rList->list[i].name);
        fprintf(stdout, "\t\tDesc : %s\n", rList->list[i].description);
        fprintf(stdout, "\t\tFlags:");
        print_beagle_flags(rList->list[i].supportFlags);
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}

void print_beagle_instance_details(BeagleInstanceDetails *inst)
{
    int rNumber = inst->resourceNumber;
    fprintf(stdout, "Using resource %i:\n", rNumber);
    fprintf(stdout, "\tRsrc Name : %s\n",inst->resourceName);
    fprintf(stdout, "\tImpl Name : %s\n", inst->implName);
    fprintf(stdout, "\tImpl Desc : %s\n", inst->implDescription);
    fprintf(stdout, "\tFlags:");
    fflush(stdout);
    print_beagle_flags(inst->flags);
    fprintf(stdout, "\n\n");
    fflush(stdout);
}

int create_beagle_instance(t_tree *tree, int quiet)
{
    if(!quiet){
        print_beagle_resource_list();
    }
    BeagleInstanceDetails inst_d;
    int num_rate_catg = tree->mod->ras->n_catg;
    int num_partials = (tree->n_otu + (tree->n_otu-2)); //taxa+internal nodes
//    int num_partials = 2*tree->n_otu-3;
    int num_scales = 2 + num_rate_catg;
    int num_branches = 2*tree->n_otu-3;
    DUMP_I(tree->n_otu, num_rate_catg, num_partials, num_branches, tree->mod->ns, tree->n_pattern, tree->mod->whichmodel);fflush(stderr);
    int beagle_inst = beagleCreateInstance(
                                  tree->n_otu,                /**< Number of tip data elements (input) */
                                  num_partials,               /**< Number of partial buffer (input) */
                                  tree->n_otu,                /**< Number of compact state representation buffers to create (input) */
                                  tree->mod->ns,              /**< Number of states in the continuous-time Markov chain (input) */
                                  tree->n_pattern,            /**< Number of site patterns to be handled by the instance (input) */
                                  1,                          /**< Number of rate matrix eigen-decomposition buffers. We supply BEAGLE the P matrix directly */
                                  num_branches,               /**< Number of rate matrix buffers (input) */
                                  num_rate_catg,              /**< Number of rate categories (input) */
                                  num_scales,                 /**< Number of scaling buffers */
                                  NULL,                       /**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
                                  0,			    /**< Length of resourceList list (input) */
                                  BEAGLE_FLAG_FRAMEWORK_CPU | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_SCALING_MANUAL,
                                  0,                /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
                                  &inst_d);
    if (beagle_inst < 0){
        fprintf(stderr, "beagleCreateInstance() failed:%i\n\n",beagle_inst);
        return beagle_inst;
    }

    if(!quiet){
        fprintf(stdout, "\nUnique BEAGLE instance id:%i\n", beagle_inst);
        print_beagle_instance_details(&inst_d);
    }

    //Set the tips
    for(int i=0; i<tree->n_otu; ++i)
    {
        assert(tree->data->c_seq[i]->len == tree->n_pattern); // number of compacts sites == number of distinct site patterns
        double* tip = tips_to_partials_nucl(tree->data->c_seq[i]->state);
        int ret = beagleSetTipPartials(beagle_inst, i, tip);
        if(ret<0){
            fprintf(stderr, "beagleSetTipPartials() on instance %i failed:%i\n\n",beagle_inst,ret);
            Free(tip);
            return ret;
        }
        Free(tip);
    }

    //Set the equilibrium freqs
    assert(tree->mod->e_frq->pi->len == tree->mod->ns);
    int ret = beagleSetStateFrequencies(beagle_inst, 0, tree->mod->e_frq->pi->v);
    if(ret<0){
        fprintf(stderr, "beagleSetStateFrequencies() on instance %i failed:%i\n\n",beagle_inst,ret);
        return ret;
    }

//    //Set the pattern weights
//    ret = beagleSetPatternWeights(beagle_inst, (const double*)tree->data->wght);
//    if(ret<0){
//        fprintf(stderr, "beagleSetPatternWeights() on instance %i failed:%i\n\n",beagle_inst,ret);
//        return ret;
//    }

//    //Set initial substitution rates (weighted)
//    ret = beagleSetCategoryRates(beagle_inst, tree->mod->ras->gamma_rr->v);
//    if(ret<0){
//        fprintf(stderr, "beagleSetCategoryRates() on instance %i failed:%i\n\n",beagle_inst,ret);
//        return ret;
//    }
//    ret = beagleSetCategoryWeights(beagle_inst, 0, /*TODO*/);
//    if(ret<0){
//        fprintf(stderr, "beagleSetCategoryWeights() on instance %i failed:%i\n\n",beagle_inst,ret);
//        return ret;
//    }

    return beagle_inst;
}

/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/
void update_beagle_partials(t_tree* tree, t_edge* b, t_node* d)
{
    /*
               |
               |<- b
               |
               d
              / \
          b1 /   \ b2
            /     \
        n_v1     n_v2
    */

    assert(!d->tax); //Partial likelihoods are only calculated on internal nodes

    //Determine d's "left" and "right" neighbors.
    t_node *n_v1, *n_v2;//d's "left" and "right" neighbor nodes
    phydbl *p_lk,*p_lk_v1,*p_lk_v2;
    phydbl *Pij1,*Pij2;
    int *sum_scale, *sum_scale_v1, *sum_scale_v2;
    int *p_lk_loc;
    n_v1 = n_v2                 = NULL;
    p_lk = p_lk_v1 = p_lk_v2    = NULL;
    Pij1 = Pij2                 = NULL;
    sum_scale_v1 = sum_scale_v2 = NULL;
    p_lk_loc                    = NULL;
    Set_All_P_Lk(&n_v1,&n_v2,
                 &p_lk,&sum_scale,&p_lk_loc,
                 &Pij1,&p_lk_v1,&sum_scale_v1,
                 &Pij2,&p_lk_v2,&sum_scale_v2,
                 d,b,tree);

    //Determine b1 and b2
    //TODO: Find a better way to do this after you understanding the index/numbering scheme for the edges
    t_edge *b1, *b2;
    b1 = b2 = NULL;
    for(int i=0;i<3;++i)//each node has 3 branches
    {
        if(NULL==b1)
            if(n_v1->b[i] == d->b[0] || n_v1->b[i] == d->b[1] || n_v1->b[i] == d->b[2])
                b1 = n_v1->b[i];
        if(NULL==b2)
            if(n_v2->b[i] == d->b[0] || n_v2->b[i] == d->b[1] || n_v2->b[i] == d->b[2])
                b2 = n_v2->b[i];
    }
//    DUMP_I(d->num, n_v1->num, n_v2->num, b->num, b1->num, b2->num);

    //Create the corresponding BEAGLE operation
    BeagleOperation operations[1] = {{d->num, BEAGLE_OP_NONE, BEAGLE_OP_NONE, n_v1->num, b1->num, n_v2->num, b2->num}};
//    BeagleOperation operations[1] = {{d->num, BEAGLE_OP_NONE, BEAGLE_OP_NONE, n_v1->tax?b1->num:n_v1->num, b1->num, n_v2->tax?b2->num:n_v2->num, b2->num}};
    int ret = beagleResetScaleFactors(tree->b_inst, 0);
    if(ret<0){
        fprintf(stderr, "beagleResetScaleFactors() on instance %i failed:%i\n\n",tree->b_inst,ret);
        Exit("");
    }
    ret = beagleUpdatePartials(tree->b_inst, operations, 1, BEAGLE_OP_NONE);
    if(ret<0){
        fprintf(stderr, "beagleUpdatePartials() on instance %i failed:%i\n\n",tree->b_inst,ret);
        Exit("");
    }

    //Fetch and Set the updated partial likelihoods
    ret = beagleGetPartials(tree->b_inst, d->num, BEAGLE_OP_NONE, (double*)p_lk);
    if(ret<0){
        fprintf(stderr, "beagleGetPartials() on instance %i failed:%i\n\n",tree->b_inst,ret);
        Exit("");
    }

    //Scaling (BEAGLE specific). That is, the p_lk vector (which is returned by BEAGLE) is indexed based on the memory layout rates * patterns * state
//    int n_patterns = tree->n_pattern;
//    phydbl p_lk_lim_inf = (phydbl)P_LK_LIM_INF;
//    int dim1 = n_patterns * tree->mod->ns;
//    int dim2 = tree->mod->ns;
//    int sum_scale_v1_val = 0;
//    int sum_scale_v2_val = 0;
//    phydbl curr_scaler;
//    phydbl smallest_p_lk;
//    int curr_scaler_pow, piecewise_scaler_pow;
//    curr_scaler = .0;
//    curr_scaler_pow = piecewise_scaler_pow = 0;
//    for(int site=0;site<n_patterns;++site)
//    {
//        for(int catg=0;catg<tree->mod->ras->n_catg;++catg)
//        {
//            smallest_p_lk  =  BIG;

//            for(int i=0;i<tree->mod->ns;++i)
//            {
//                if(p_lk[catg*dim1+site*dim2+i] < smallest_p_lk)
//                    smallest_p_lk = p_lk[catg*dim1+site*dim2+i];
//            }

//            sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
//            sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);

//            sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
////            DUMP_D(smallest_p_lk, curr_scaler);
////            DUMP_I(curr_scaler_pow);
////            fprintf(stderr,"\n%d",sum_scale[catg*n_patterns+site]);

//            /* Scaling */
//            if(smallest_p_lk < p_lk_lim_inf)
//              {
//                curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
//                curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);

//                sum_scale[catg*n_patterns+site] += curr_scaler_pow;

//                do
//                  {
//                    piecewise_scaler_pow = MIN(curr_scaler_pow,63);
//                    curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
//                    for(int j=0;j<tree->mod->ns;++j)
//                      {
//                        p_lk[catg*dim1+site*dim2+j] *= curr_scaler;

//                        if(p_lk[catg*dim1+site*dim2+j] > BIG)
//                          {
//                            PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
//                            PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
//                            PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
//                            Warn_And_Exit("\n");
//                          }
//                      }
//                    curr_scaler_pow -= piecewise_scaler_pow;
//                  }
//                while(curr_scaler_pow != 0);
//              }
//        }
//    }

//    int site1,catg1;
//    for(site1=0;site1<n_patterns;++site1)
//    {
//        for(catg1=0;catg1<tree->mod->ras->n_catg;++catg1)
//        {
//            int i;
//            for(i=0;i<tree->mod->ns;++i)
//            {
//                fprintf(stdout, "%f\n",p_lk[catg1*dim1+site1*dim2+i]);
//            }
//        }
//    }

//    int parent_nodes[2] = {d->num, d->num};//The child nodes share the same parent
//    int child_nodes[2] = {n_v1->num, n_v2->num};
//    int p_mat[2] = {b1->num, b2->num};
//    int category_weights = 0;
//    int state_freqs = 0;
//    int cum_scale = BEAGLE_OP_NONE;
//    double* log_lks = (double*)malloc(2*sizeof(double)); if(NULL==log_lks) Warn_And_Exit(__PRETTY_FUNCTION__);
//    ret = beagleCalculateEdgeLogLikelihoods(tree->b_inst, &parent_nodes, &child_nodes, &p_mat, NULL, NULL, &category_weights, &state_freqs, &cum_scale, 2, log_lks, NULL, NULL);
//    if(ret<0){
//        fprintf(stderr, "beagleCalculateEdgeLogLikelihoods() on instance %i failed:%i\n\n",tree->b_inst,ret);
//        Exit("");
//    }
}

int finalize_beagle_instance(t_tree *tree)
{
    if(tree->b_inst > 0)
    {
        int ret = beagleFinalizeInstance(tree->b_inst);
        if(ret<0) fprintf(stderr, "\nFailed to finalize BEAGLE instance %i: %i\n\n", tree->b_inst, ret);
        return ret;
    }
    return 0;
}


#endif // BEAGLE_UTILS_CPP

