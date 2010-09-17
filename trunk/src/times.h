/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef TIMES_H
#define TIMES_H

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"
#include "mcmc.h"
#include "stats.h"
#include "tiporder.h"
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

int  TIMES_main(int argc, char **argv);
void TIMES_Bl_From_T_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void TIMES_Bl_From_T(t_tree *tree);
void TIMES_Optimize_Node_Times_Serie(t_node *a, t_node *d, t_tree *tree);
void TIMES_Round_Optimize(t_tree *tree);
void TIMES_Print_Node_Times(t_node *a, t_node *d, t_tree *tree);
t_edge *TIMES_Find_Best_Root_Position(t_tree *tree);
void TIMES_Least_Square_Node_Times(t_edge *e_root, t_tree *tree);
void TIMES_Least_Square_Node_Times_Pre(t_node *a, t_node *d, phydbl *A, phydbl *b, int n, t_tree *tree);
void TIMES_Mult_Time_Stamps(t_tree *tree);
void TIMES_Div_Time_Stamps(t_tree *tree);
void TIMES_Optimize_Tree_Height(t_tree *tree);
void TIMES_Adjust_Node_Times(t_tree *tree);
void TIMES_Adjust_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree);
void TIMES_Optimize_Root_Height(t_tree *tree);
void TIMES_Estimate_Branch_Rates(t_tree *tree);
t_edge *TIMES_Find_Best_Root_Position_Approx(t_tree *tree);
void TIMES_Estimate_Branch_Rate_Parameter(t_tree *tree);
phydbl TIMES_Classify_Branch_In_Rate_Class(t_tree *tree);
void TIMES_Compute_Rates_And_Times_Least_Square_Adjustments(t_tree *tree);
void TIMES_Compute_Rates_And_Times_Least_Square_Adjustments_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree);
void TIMES_Classify_Branch_Rates(t_tree *tree);
int TIMES_Check_MC(t_tree *tree);
void TIMES_Set_All_Node_Priors(t_tree *tree);
void TIMES_Set_All_Node_Priors_Bottom_Up(t_node *a, t_node *d, t_tree *tree);
void TIMES_Set_All_Node_Priors_Top_Down(t_node *a, t_node *d, t_tree *tree);

#endif
