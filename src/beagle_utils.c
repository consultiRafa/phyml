#ifndef BEAGLE_UTILS_CPP
#define BEAGLE_UTILS_CPP

#include  <stdio.h>
#include "beagle_utils.h"


void printFlags(long inFlags) {
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
        printFlags(rList->list[i].supportFlags);
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}

void print_beagle_instance_details(BeagleInstanceDetails *inst)
{
    int rNumber = inst->resourceNumber;
    fprintf(stdout, "Using resource %i:\n", rNumber);
    fprintf(stdout, "\tRsrc Name : %s\n",inst->resourceName);
    fprintf(stdout, "\tImpl Name : %s\n", inst->implName);
    fprintf(stdout, "\tImpl Desc : %s\n", inst->implDescription);
    fprintf(stdout, "\tFlags:");
    printFlags(inst->flags);
    fprintf(stdout, "\n\n");
}

int create_beagle_instance(t_tree *tree, int quiet)
{
    if(!quiet)  print_beagle_resource_list();
    BeagleInstanceDetails inst_d;
    int num_rate_catg = tree->mod->ras->n_catg;
    int num_partials = (tree->n_otu-2) + num_rate_catg;  //Number of internal nodes + Number of roots
    int num_scales = 2 + num_rate_catg;
    int num_branches = 2*tree->n_otu-3;
    DUMP_I(tree->n_otu, num_rate_catg, num_partials, num_scales, num_branches, tree->mod->ns, tree->n_pattern );
    int beagle_inst = beagleCreateInstance(
                                  tree->n_otu,                /**< Number of tip data elements (input) */
                                  num_partials,               /**< Number of partial buffer (input) */
                                  tree->n_otu,                /**< Number of compact state representation buffers to create (input) */
                                  tree->mod->ns,              /**< Number of states in the continuous-time Markov chain (input) */
                                  tree->n_pattern,            /**< Number of site patterns to be handled by the instance (input) */
                                  1,                          /**< Number of rate matrix eigen-decomposition buffers to allocate (input) */
                                  num_branches,               /**< Number of rate matrix buffers (input) */
                                  num_rate_catg,              /**< Number of rate categories (input) */
                                  num_scales,                 /**< Number of scaling buffers */
                                  NULL,			/**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
                                  0,			    /**< Length of resourceList list (input) */
                                  BEAGLE_FLAG_FRAMEWORK_CPU | BEAGLE_FLAG_PRECISION_SINGLE | BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_SCALING_MANUAL,
                                  0,                /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
                                  &inst_d);
    if(!quiet) fprintf(stdout, "\nUnique BEAGLE instance id:%i\n", beagle_inst);
    if (beagle_inst < 0)
    {
        fprintf(stderr, "Failed to obtain BEAGLE instance\n\n");
        return beagle_inst;
    }
    else if(!quiet) print_beagle_instance_details(&inst_d);
    return beagle_inst;
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

