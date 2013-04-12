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
#include "mixt.h"
#include "sergeii.h"
#ifdef MPI
#include "mpi_boot.h"
#endif
#include <stdlib.h>


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Function checks if the randomized node times are within the 
//upper and lower time limits, taken into account the times of 
//the ancestor and descendent.
void Check_Node_Time(t_node *a, t_node *d, int *result, t_tree *tree)
{
  phydbl t_low, t_up;
  phydbl *t_prior_min, *t_prior_max, *nd_t;

  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;

  if(a == tree -> n_root && (nd_t[a -> num] > MIN(t_prior_max[a -> num], MIN(nd_t[a -> v[1] -> num], nd_t[a -> v[2] -> num]))))  
    {
      *result = FALSE;
      return;
    }
  if(d -> tax) return;
  else
    { 
      t_low = MAX(t_prior_min[d -> num], nd_t[d -> anc -> num]);
      t_up = MIN(t_prior_max[d -> num], MIN(nd_t[d -> v[1] -> num], nd_t[d -> v[2] -> num]));
      if(nd_t[d -> num] < t_low || nd_t[d -> num] > t_up)
        {
          *result = FALSE; 
        }

      int i;
      For(i,3) 
	if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root))
             Check_Node_Time(d, d -> v[i], result, tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Function calculates the TOTAL number of calibration combinations, 
//given the number of nodes to which each calibartion applies to.
int Number_Of_Comb(t_cal *calib)
{

  int num_comb;

  if(!calib) return(1);
  num_comb = 1;
  do
    {
      num_comb *= calib -> n_all_applies_to;
      if(calib -> next) calib = calib -> next;
      else break;
    }
  while(calib);
  return(num_comb);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Function sets current calibartion in the following way:
//Suppose we have a vector of calibrations C=(C1, C2, C3), each calibration  
//applies to a set of nodes. we can reach each node number through the indeces (corresponds 
//to the number the information was read). C1={0,1,2}, C2={0,1}, C3={0};
//The total number of combinations is 3*2*1=6. The first combination with row number 0 
//will be {0,0,0}, the second row will be {0,1,0} and so on. Calling the node numbers with 
//the above indeces will return current calibration.   
void Set_Current_Calibration(int row, t_tree *tree)
{

  t_cal *calib;
  phydbl *t_prior_min, *t_prior_max; 
  short int *t_has_prior;
  int k, i, j, *curr_nd_for_cal;

  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  curr_nd_for_cal = tree -> rates -> curr_nd_for_cal;

  for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
    {
      t_prior_min[j] = -BIG;
      t_prior_max[j] = BIG;
      t_has_prior[j] = NO; 
    }
  
  k = -1;
  i = 0;
  do
    {
      k = (row % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);      
      t_prior_min[calib -> all_applies_to[k] -> num] = calib -> lower;
      t_prior_max[calib -> all_applies_to[k] -> num] = calib -> upper;
      t_has_prior[calib -> all_applies_to[k] -> num] = YES; 
      curr_nd_for_cal[i] = calib -> all_applies_to[k] -> num;
      i++;
      if(calib->next) calib = calib->next;
      else break;    
    }
  while(calib);
  //while(calib -> prev) calib = calib -> prev;
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Calculate the prior probability for node times taking into account the 
//probailitis with which each calibration applies to the particular node.
phydbl TIMES_Calib_Cond_Prob(t_tree *tree)
{

  phydbl times_lk, *Yule_val, *times_partial_proba, times_tot_proba, *t_prior_min, *t_prior_max, min_value, c, constant, ln_t; 
  short int *t_has_prior;
  int i, j, k, tot_num_comb, result;
  t_cal *calib;
 

  times_tot_proba = 0.0;
  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  //constant = tree -> K;
  tot_num_comb = Number_Of_Comb(calib);

  
  Yule_val = (phydbl *)mCalloc(tot_num_comb, sizeof(phydbl));    
  times_partial_proba = (phydbl *)mCalloc(tot_num_comb, sizeof(phydbl));   

  For(i, tot_num_comb)
    {
      for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
        {
          t_prior_min[j] = -BIG;
          t_prior_max[j] = BIG;
          t_has_prior[j] = NO; 
        }
      times_partial_proba[i] = 1.0; 
      do
        {
          k = (i % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);
          t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower);
          t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper);
          t_has_prior[calib -> all_applies_to[k] -> num] = YES;
          if((t_prior_min[calib -> all_applies_to[k] -> num] > t_prior_max[calib -> all_applies_to[k] -> num])) times_partial_proba[i] = 0.0; 
          else times_partial_proba[i] *= calib -> proba[calib -> all_applies_to[k] -> num]; 
          if(calib -> next) calib = calib -> next;
          else break;
        }
      while(calib);
      TIMES_Set_All_Node_Priors(tree);
      //printf("\n. p[%i] = %f \n", i + 1, times_partial_proba[i]);
      result = TRUE;
      Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree) ;
      Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree) ;
      if(result != TRUE) times_partial_proba[i] = 0.0; 
      
      //tree -> rates -> birth_rate = 4.0;
      
      times_lk = TIMES_Lk_Yule_Order(tree);
      constant = 1.0; 
      
      if(times_lk > -INFINITY) constant = Slicing_Calibrations(tree);

      //printf("\n. K = [%f] \n", K);
      //K = Norm_Constant_Prior_Times(tree);
      //printf("\n. K = [%f] \n", K[0]); printf("\n. K = [%f] \n", K[1]);
      //Exit("\n");
      //Yule_val[i] = K[i] * TIMES_Lk_Yule_Order(tree);
      //For(j, 2 * tree -> n_otu - 1) printf("\n. [2] Node [%d] min [%f] max[%f]\n", j, tree -> rates -> t_prior_min[j], tree -> rates -> t_prior_max[j]);

      Yule_val[i] = constant * times_lk;

      //Yule_val[i] = TIMES_Lk_Yule_Order(tree);
 
      while(calib -> prev) calib = calib -> prev;
    }
 
  min_value = 0.0;
  For(i, tot_num_comb) if(Yule_val[i] < min_value && Yule_val[i] > -INFINITY) min_value = Yule_val[i];

  c = -600. - min_value;
  times_tot_proba = 0.0;
  For(i, tot_num_comb)
    {
      times_tot_proba += times_partial_proba[i] * EXP(Yule_val[i] + c);
    }

  For(i, 2 * tree -> n_otu - 1) t_has_prior[i] = NO;

  ln_t = -c + LOG(times_tot_proba);

  free(Yule_val);  
  free(times_partial_proba);

  return(ln_t);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Randomly choose a combination of calibrations drawing an index of calibration combination, 
//used function Set_Cur_Calibration.
void Random_Calibration(t_tree *tree)
{
  int rnd, num_comb;
  t_cal *calib;

  calib = tree -> rates -> calib;

  num_comb =  Number_Of_Comb(calib);

  srand(time(NULL));
  rnd = rand()%(num_comb);
  
  Set_Current_Calibration(rnd, tree);
  TIMES_Set_All_Node_Priors(tree); 

}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Variable curr_nd_for_cal is a vector of node numbers, the length of that vector is a number of calibrations. 
//Function randomly updates that vector by randomly changing one node and setting times limits with respect 
//to a new vector.  
int RND_Calibration_And_Node_Number(t_tree *tree)
{
  int i, j, tot_num_cal, cal_num, node_ind, node_num, *curr_nd_for_cal;
  phydbl *t_prior_min, *t_prior_max;
  short int *t_has_prior;
  t_cal *cal;

  tot_num_cal = tree -> rates -> tot_num_cal;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  curr_nd_for_cal = tree -> rates -> curr_nd_for_cal;
  cal = tree -> rates -> calib;

  cal_num = rand()%(tot_num_cal - 1);
    
  i = 0;
  while (i != cal_num)
    {
      cal = cal -> next;
      i++;
    }

  node_ind = rand()%(cal -> n_all_applies_to);
  node_num = cal -> all_applies_to[node_ind] -> num;

  curr_nd_for_cal[cal_num] = node_num;

  for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
    {
      t_prior_min[j] = BIG;
      t_prior_max[j] = BIG;
      t_has_prior[j] = NO; 
    }

  while(cal -> prev) cal = cal -> prev;
  
  i = 0;
  do
    {
      t_prior_min[curr_nd_for_cal[i]] = cal -> lower;
      t_prior_max[curr_nd_for_cal[i]] = cal -> upper;
      t_has_prior[curr_nd_for_cal[i]] = YES; 
      i++;
      if(cal->next) cal = cal -> next;
      else break;    
    }
  while(cal);

  while(cal -> prev) cal = cal -> prev;

  TIMES_Set_All_Node_Priors(tree);

  return(node_num);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Return the value uniformly distributed between two values.
phydbl Randomize_One_Node_Time(phydbl min, phydbl max)
{
  phydbl u;

  u = Uni();
  u *= (max - min);
  u += min;

  return(u);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Calculates the Hastings ratio for the analysis. Used in case of 
//calibration conditional jump. NOT THE RIGHT ONE TO USE!
void Lk_Hastings_Ratio_Times(t_node *a, t_node *d, phydbl *tot_prob, t_tree *tree)
{
  phydbl t_low, t_up;
  phydbl *t_prior_min, *t_prior_max, *nd_t;

  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;

  if(d -> tax) return;
  else
    { 
      t_low = MAX(t_prior_min[d -> num], nd_t[d -> anc -> num]);
      t_up = MIN(t_prior_max[d -> num], MIN(nd_t[d -> v[1] -> num], nd_t[d -> v[2] -> num]));

      (*tot_prob) += LOG(1) - LOG(t_up - t_low);

      int i;
      For(i,3) 
	if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root))
          { 
              Lk_Hastings_Ratio_Times(d, d -> v[i], tot_prob, tree);
          }
    } 
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Updates nodes which are below a randomized node in case if new proposed time 
//for that node is below the current value. 
void Update_Descendent_Cond_Jump(t_node *a, t_node *d, phydbl *L_Hast_ratio, t_tree *tree)
{
  int result = TRUE;
  phydbl t_low, t_up;
  phydbl *t_prior_min, *t_prior_max, *nd_t;

  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;

  Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree);  
  Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree);

  if(d -> tax) return;
  else
  {
    if(result != TRUE)
      {  
        int i;
        t_low = MAX(nd_t[a -> num], t_prior_min[d -> num]);
        if(t_low < MIN(nd_t[d -> v[1] -> num], nd_t[d -> v[2] -> num])) t_up  = MIN(t_prior_max[d -> num], MIN(nd_t[d -> v[1] -> num], nd_t[d -> v[2] -> num])); 
        else t_up  = t_prior_max[d -> num]; 
        nd_t[d -> num] = Randomize_One_Node_Time(t_low, t_up);
        (*L_Hast_ratio) += LOG(1) - LOG(t_up - t_low);
        For(i,3) 
          if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root)) 
            Update_Descendent_Cond_Jump(d, d -> v[i], L_Hast_ratio, tree);
      }
    else return;
  }
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Updates nodes which are above a randomized node in case if new proposed time 
//for that node is above the current value.
void Update_Ancestor_Cond_Jump(t_node *d, phydbl *L_Hast_ratio, t_tree *tree)
{
  int result = TRUE;
  phydbl t_low, t_up;
  phydbl *t_prior_min, *t_prior_max, *nd_t;

  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;
  
  Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree);  
  Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree);

  if(result != TRUE) 
    {      
      if(d == tree -> n_root)
        {
                      
          t_low = t_prior_min[d -> num];
          t_up  = MIN(t_prior_max[d -> num], MIN(nd_t[d -> v[1] -> num], nd_t[d -> v[2] -> num])); 
          nd_t[d -> num] = Randomize_One_Node_Time(t_low, t_up);
          (*L_Hast_ratio) += LOG(1) - LOG(t_up - t_low);
          return;
        }
      else
        { 
          t_up  = MIN(t_prior_max[d -> num], MIN(nd_t[d -> v[1] -> num], nd_t[d -> v[2] -> num]));           
          if(nd_t[d -> anc -> num] > t_up) t_low =  t_prior_min[d -> num];
          else  t_low  = MAX(t_prior_min[d -> num], nd_t[d -> anc -> num]); 
          nd_t[d -> num] = Randomize_One_Node_Time(t_low, t_up);
          (*L_Hast_ratio) += LOG(1) - LOG(t_up - t_low);
          Update_Ancestor_Cond_Jump(d -> anc, L_Hast_ratio, tree); 
        }
      
    }
  else return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//when made a calibration conditional jump, updates node times 
//with respect to the new calibration which was made with respect 
//to the randomly chosen node, the root is fixed. Updates only those nodes 
//that are not within new intervals. Traverse up and down.
void Update_Times_RND_Node_Ancestor_Descendant(int rnd_node, phydbl *L_Hast_ratio, t_tree *tree)
{
  int i;
  phydbl *t_prior_min, *t_prior_max, *nd_t;
  phydbl new_time_rnd_node = 0.0;

  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;

  new_time_rnd_node = Randomize_One_Node_Time(t_prior_min[rnd_node], t_prior_max[rnd_node]);

  nd_t[rnd_node] = new_time_rnd_node; 

  Update_Ancestor_Cond_Jump(tree -> a_nodes[rnd_node] -> anc, L_Hast_ratio, tree);
  For(i,3) 
    if((tree -> a_nodes[rnd_node] -> v[i] != tree -> a_nodes[rnd_node] -> anc) && (tree -> a_nodes[rnd_node] -> b[i] != tree -> e_root)) 
      Update_Descendent_Cond_Jump(tree -> a_nodes[rnd_node], tree -> a_nodes[rnd_node] -> v[i], L_Hast_ratio, tree);
 
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//when made a calibration conditional jump, updates node times 
//with respect to the new calibration which was made with respect 
//to the randomly chosen node, starting from the root down to the tips. 
//Updates only those nodes that are not within new intervals. 
void Update_Times_Down_Tree(t_node *a, t_node *d, phydbl *L_Hastings_ratio, t_tree *tree)
{
  int i; 
  phydbl *t_prior_min, *t_prior_max, *nd_t, t_low, t_up;

  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;


  t_low = MAX(t_prior_min[d -> num], nd_t[d -> anc -> num]);
  t_up  = t_prior_max[d -> num];  

  d -> anc = a;
  //printf("\n. [1] Node number: [%d] \n", d -> num);
  if(d -> tax) return;
  else
  {    
    if(nd_t[d -> num] > t_up || nd_t[d -> num] < t_low)
      { 
        //printf("\n. [2] Node number: [%d] \n", d -> num);
        //(*L_Hastings_ratio) += (LOG(1) - LOG(t_up - t_low));
        (*L_Hastings_ratio) += (- LOG(t_up - t_low));
        nd_t[d -> num] = Randomize_One_Node_Time(t_low, t_up);
        t_prior_min[d -> num] = t_low;
        t_prior_max[d -> num] = t_up;
      }    
   
    For(i,3) 
      if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root)) 
        Update_Times_Down_Tree(d, d -> v[i], L_Hastings_ratio, tree);
  }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void PhyTime_XML(char *xml_file)
{

  FILE *f;
  char **clade, *clade_name, **mon_list;
  phydbl low, up;
  int i, j, n_taxa, clade_size, node_num, n_mon; //rnd_num
  xml_node *n_r, *n_t, *n_m, *n_cur;
  t_cal *last_calib; 
  /* t_cal *cur; */
  align **data, **c_seq;
  option *io;
  calign *cdata;
  t_opt *s_opt;
  t_mod *mod;
  time_t t_beg,t_end;
  int r_seed;
  char *most_likely_tree;
  int user_lk_approx;
  t_tree *tree;
  t_node **a_nodes; //*node;
  m4 *m4mod;
 
  srand(time(NULL));

  i = 0;
  j = 0;
  last_calib       = NULL;
  mod              = NULL;
  most_likely_tree = NULL;
  n_taxa = 0; 

  //file can/cannot be open:
  if ((f =(FILE *)fopen(xml_file, "r")) == NULL)
    {
      PhyML_Printf("\n==File can not be open...\n");
      Exit("\n");
    }

  n_r = XML_Load_File(f);

  //memory allocation for model parameters:
  io = (option *)Make_Input();
  mod   = (t_mod *)Make_Model_Basic();
  s_opt = (t_opt *)Make_Optimiz();
  m4mod = (m4 *)M4_Make_Light();
  Set_Defaults_Input(io);                                                                          
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);  
  io -> mod = mod;
  mod = io -> mod;
  mod -> s_opt = s_opt;
  clade_size = -1;

  ////////////////////////////////////////////////////////////////////////////
  //////////////////////reading tree topology:////////////////////////////////

  //looking for a node <topology>
  n_t = XML_Search_Node_Name("topology", YES, n_r);

  //setting tree:
  tree        = (t_tree *)mCalloc(1,sizeof(t_tree));
  n_cur = XML_Search_Node_Name("instance", YES, n_t);
  if(n_cur != NULL)
    {
      if(XML_Search_Attribute(n_cur, "user.tree") != NULL) 
	{
	  strcpy(io -> out_tree_file, XML_Search_Attribute(n_cur, "user.tree") -> value); 
	  io -> fp_out_tree  = Openfile(io -> out_tree_file, 1);
	  io -> tree = Read_Tree_File_Phylip(io -> fp_in_tree);
	}
      else
	{
	  PhyML_Printf("\n==Tree was not found. \n");
	  PhyML_Printf("\n==Either specify tree file name or enter the whole tree. \n");
	  Exit("\n");
	}
    }
  else io -> tree  = Read_Tree(&n_t -> value);
  io -> n_otu = io -> tree -> n_otu;
  tree = io -> tree;


  //setting initial values to n_calib:
  For(i, 2 * tree -> n_otu - 2) 
    {
      tree -> a_nodes[i] -> n_calib = 0;
      //PhyML_Printf("\n. '%d' \n", tree -> a_nodes[i] -> n_calib);
    }


  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //memory for nodes:
  a_nodes   = (t_node **)mCalloc(2 * io -> n_otu - 1,sizeof(t_node *));
  For(i, 2 * io -> n_otu - 2) a_nodes[i] = (t_node *)mCalloc(1,sizeof(t_node));
 
  //setting a model: 
  tree -> rates = RATES_Make_Rate_Struct(io -> n_otu);                                    
  RATES_Init_Rate_Struct(tree -> rates, io -> rates, tree -> n_otu);

  //reading seed:
  if(XML_Search_Attribute(n_r, "seed")) io -> r_seed = String_To_Dbl(XML_Search_Attribute(n_r, "seed") -> value);
   
  //TO DO: check that the tree has a root...
  Update_Ancestors(io -> tree -> n_root, io -> tree -> n_root -> v[2], io -> tree);
  Update_Ancestors(io -> tree -> n_root, io -> tree -> n_root -> v[1], io -> tree);
		  

  ////////////////////////////////////////////////////////////////////////////
  //////////////////////memory allocation for temp parameters/////////////////

  //memory for monitor flag:
  io -> mcmc -> monitor = (int *)mCalloc(2 * io -> n_otu - 1,sizeof(int));

 
  //memory for sequences:
  n_cur = XML_Search_Node_Name("alignment", YES, n_r);

  data   = (align **)mCalloc(io -> n_otu,sizeof(align *));
  For(i, io -> n_otu)
    {
      data[i]          = (align *)mCalloc(1,sizeof(align));
      data[i] -> name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      if(n_cur -> child -> value != NULL) data[i] -> state = (char *)mCalloc(strlen(n_cur -> child -> value) + 1,sizeof(char));
      else data[i] -> state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
    }
  io -> data = data;
  //tree -> data = data;

 
  //memory for clade:
  clade_name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  clade = (char **)mCalloc(T_MAX_FILE,sizeof(char *));
  For(i, T_MAX_FILE)
    {
      clade[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  //memory for list of clades to be monitored
  mon_list = (char **)mCalloc(T_MAX_FILE,sizeof(char *));
  For(i, T_MAX_FILE)
    {
      mon_list[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 
  //reading monitor node:
  i = 0;
  n_m = XML_Search_Node_Name("monitor", YES, n_r);
  if(n_m != NULL)
    {
      do
	{
	  strcpy(mon_list[i], n_m -> child -> attr -> value);
	  i++;
	  if(n_m -> child) n_m -> child = n_m -> child -> next;
	  else break;
	}
      while(n_m -> child);
      n_mon = i;
    }
  else
    {
      n_mon = 0;
      PhyML_Printf("\n. There is no clade to be monitored. \n");
    }

 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////

  //chekcing for calibration node (upper or lower bound) to exist:
  n_cur = XML_Search_Node_Name("calibration", YES, n_r);

  if(n_cur == NULL)
    {
      PhyML_Printf("\n==There is no calibration information provided. \n");
      PhyML_Printf("\n==Please check your data. \n");
      Exit("\n");
    }
  else
    {
      if(XML_Search_Node_Name("upper", NO, n_cur -> child) == NULL && XML_Search_Node_Name("lower", NO, n_cur -> child) == NULL)
	{
	  PhyML_Printf("\n==There is no calibration information provided. \n");
	  PhyML_Printf("\n==Please check your data. \n");
	  Exit("\n");
	}
    }

 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
  if(XML_Search_Attribute(n_r, "use_data") != NULL) 
    {
      if(XML_Search_Attribute(n_r, "use_data") -> value != NULL && (!strcmp(XML_Search_Attribute(n_r, "use_data") -> value, "YES")))
        io -> mcmc -> use_data = YES;
      if(XML_Search_Attribute(n_r, "use_data") -> value != NULL && (!strcmp(XML_Search_Attribute(n_r, "use_data") -> value, "NO")))
        io -> mcmc -> use_data = NO;
    }

  n_r = n_r -> child; 
  tree -> rates -> tot_num_cal = 0;
  do
    {
      if(!strcmp(n_r -> name, "alignment"))//looking for a node <alignment>.
	{
	  if(!n_r -> attr -> value)
	  {
		PhyML_Printf("\n==Not found sequence type (nt / aa). \n");
		PhyML_Printf("\n==Please, include data to node <%s> attribute value. \n", n_r -> name);
		Exit("\n");
	  }	
	  if(!strcmp(To_Upper_String(n_r -> attr -> value), "NT")) 
	    {
	      io -> datatype = 0;
	      io -> mod -> ns = 4;
	    }
 	  if(!strcmp(To_Upper_String(n_r -> attr -> value), "AA")) 
	    {
	      io -> datatype = 1;
	      io -> mod -> ns = 20;
	    }

	  n_cur = XML_Search_Node_Name("instance", YES, n_r);
	  if(n_cur != NULL)
	    {
	      if(XML_Search_Attribute(n_cur, "user.alignment") != NULL) 
		{
		  strcpy(io -> in_align_file, XML_Search_Attribute(n_cur, "user.alignment") -> value); 
		  io -> fp_in_align  = Openfile(io -> in_align_file, 1);
		  Detect_Align_File_Format(io);
		  io -> data = Get_Seq(io);
		}
	    }
	  else
	    {	
	      i = 0;
	      do
		{
                  strcpy(io -> in_align_file, "sergeii"); 
		  strcpy(io -> data[i] -> name, n_r -> child -> attr -> value);
		  strcpy(io -> data[i] -> state, To_Upper_String(n_r -> child -> value));
		  i++;
		  if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
		  else break;
		}
	      while(n_r -> child);
	      n_taxa = i;
	    }

	  //checking if a sequences of the same lengths:

	  i = 1;
	  For(i, n_taxa) if(strlen(io -> data[0] -> state) != strlen(io -> data[i] -> state))
	    {
	      printf("\n. Sequences are of different length. Please check your data...\n");
	      Exit("\n");
	      break;
	    }  

	  //checking sequence names:
	  i = 0;
	  For(i, n_taxa) Check_Sequence_Name(io -> data[i] -> name);

	  //check if a number of tips is equal to a number of taxa:
	  if(n_taxa != io -> n_otu)
	    {
	      PhyML_Printf("\n==Number of taxa is not the same as a number of tips. Check your data...\n");
	      Exit("\n");
	    }

          //deleting '-', etc. from sequences:
          io -> data[0] -> len = strlen(io -> data[0] -> state);
          Post_Process_Data(io);   
 	  n_r = n_r -> next;
	}
      else if(!strcmp(n_r -> name, "calibration"))//looking for a node <calibration>.
	{
          tree -> rates -> tot_num_cal++;
	  if (tree -> rates -> calib == NULL) tree -> rates -> calib = Make_Calib(tree -> n_otu);
          if(last_calib)
            {
              last_calib -> next = tree -> rates -> calib;
              tree -> rates -> calib -> prev = last_calib;
            }
          last_calib = tree -> rates -> calib;

	  low = -BIG;
	  up  = BIG;
	  n_cur = XML_Search_Node_Name("lower", YES, n_r);
	  if(n_cur != NULL) low = String_To_Dbl(n_cur -> value);
	  n_cur = XML_Search_Node_Name("upper", YES, n_r);
	  if(n_cur != NULL) up = String_To_Dbl(n_cur -> value);
	  do
	    {
	      if(!strcmp("appliesto", n_r -> child -> name)) 
		{
		  strcpy(clade_name, n_r -> child -> child -> attr -> value);//reached clade names
		  if(!strcmp("@root@", clade_name))
		    {
		      node_num = io -> tree -> n_root -> num;
		      For(j, n_mon)
			{
			  if(!strcmp("@root@", mon_list[j])) io -> mcmc -> monitor[node_num] = YES;
			}	
		    }
		  else		  
		    {
		      i = 0; 
		      do
			{
			  strcpy(clade[i], n_r -> child -> child -> child -> attr -> value); 
			  i++;
			  if(n_r -> child -> child -> child -> next) n_r -> child -> child -> child = n_r -> child -> child -> child -> next;
			  else break;
			}
		      while(n_r -> child -> child -> child);
		      clade_size = i;
		      node_num = Find_Clade(clade, clade_size, io -> tree);
			{
			  if(!strcmp(clade_name, mon_list[j])) io -> mcmc -> monitor[node_num] = YES;
			}	
                    }     
                  ////////////////////////////////////////////////////////////////////////////////////////////////////
                  if(n_r -> child -> attr) tree -> rates -> calib -> proba[node_num] = String_To_Dbl(n_r -> child -> attr -> value);
                  if(!n_r -> child -> attr && n_r -> child -> next == NULL)tree -> rates -> calib -> proba[node_num] = 1.;
                  if(!n_r -> child -> attr && n_r -> child -> next)
                    {
                      PhyML_Printf("==You either need to provide information about probability with which calibration \n");
                      PhyML_Printf("==applies to a node or you need to apply calibartion only to one node. \n");
                      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
                      Exit("\n");
                    }
                  tree -> rates -> calib -> all_applies_to[tree -> rates -> calib -> n_all_applies_to] -> num = node_num; 
                  tree -> rates -> calib -> n_all_applies_to++;
                  tree -> rates -> calib -> lower = low;
                  tree -> rates -> calib -> upper = up;
                  ////////////////////////////////////////////////////////////////////////////////////////////////////
                  /*
                    if(tree -> a_nodes[node_num] -> calib == NULL) 
                    {
                      tree -> a_nodes[node_num] -> calib = (t_cal **)mCalloc(2 * tree -> n_otu - 1, sizeof(t_cal *));
                      For(i, 2 * tree -> n_otu - 1) tree -> a_nodes[node_num] -> calib[i] = Make_Calib(tree -> n_otu);
                      tree -> a_nodes[node_num] -> calib_applies_to = (short int *)mCalloc(2 * tree -> n_otu - 1, sizeof(short int));
                    }
                  tree -> a_nodes[node_num] -> calib[tree -> a_nodes[node_num] -> n_calib] -> lower = low;
                  tree -> a_nodes[node_num] -> calib[tree -> a_nodes[node_num] -> n_calib] -> upper = up;
                  tree -> a_nodes[node_num] -> calib[tree -> a_nodes[node_num] -> n_calib] -> calib_proba =  String_To_Dbl(n_r -> child -> attr -> value);
                  tree -> a_nodes[node_num] -> n_calib++;
                  */
                  /////////////////////////////////////////////////////////////////////////////////////////////////////
                  PhyML_Printf("\n. .......................................................................");
                  PhyML_Printf("\n");
                  PhyML_Printf("\n. Node number to which calibration applies to is: [%d]",node_num);
                  PhyML_Printf("\n. Clade name: [%s]", clade_name);
                  if(strcmp(clade_name, "@root@"))
                    {
                      For(i, clade_size) PhyML_Printf("\n. Taxon name: [%s]", clade[i]);
                    }
                  PhyML_Printf("\n. Lower bound set to: %15f time units.", low);
                  PhyML_Printf("\n. Upper bound set to: %15f time units.", up);
                  PhyML_Printf("\n. .......................................................................");            
                  //////////////////////////////////////////////////////////////////////////////////////////////////////
                  if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
                  else break;
		}
	      else if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
	      else break;	      
	    }
	  while(n_r -> child);                  
          //PhyML_Printf("\n. '%d'\n", tree -> rates -> calib -> n_all_applies_to);
          tree -> rates -> calib = tree -> rates -> calib -> next;	   
	  n_r = n_r -> next;
	}
      else if(!strcmp(n_r -> name, "ratematrices"))//initializing rate matrix:
        {
          if(n_r -> child) 
            {
              Make_Ratematrice_From_XML_Node(n_r -> child, io, mod);
              n_r = n_r -> next;
            } 
          else n_r = n_r -> next;
        }
      else if(!strcmp(n_r -> name, "equfreqs"))//initializing frequencies:
        {
           if(n_r -> child) 
             {
               Make_Efrq_From_XML_Node(n_r -> child , io, mod);
               n_r = n_r -> next;
             }
           else n_r = n_r -> next;
        }
      else if(!strcmp(n_r -> name, "siterates"))//initializing site rates:
        {
          if(n_r -> child) 
            {
              Make_RAS_From_XML_Node(n_r, io -> mod);
              n_r = n_r -> next;
            }
          else n_r = n_r -> next;
        }
      else if (n_r -> next) n_r = n_r -> next;
      else break;
    }
  while(1);
 //Root to be the last one calibration read from the file.. 
  n_cur = XML_Search_Node_Name("calibration_root", NO, n_r -> parent);
  if(n_cur)
    {         
      tree -> rates -> tot_num_cal++;
      if (tree -> rates -> calib == NULL) tree -> rates -> calib = Make_Calib(tree -> n_otu);
      if(last_calib)
        {
          last_calib -> next = tree -> rates -> calib;
          tree -> rates -> calib -> prev = last_calib;
        }
      last_calib = tree -> rates -> calib;
      
      low = -BIG;
      up  =  BIG;     
      n_cur = XML_Search_Node_Name("lower", YES, n_cur);
      if(n_cur != NULL) low = String_To_Dbl(n_cur -> value);
      n_cur = XML_Search_Node_Name("upper", NO, n_cur -> parent);
      if(n_cur != NULL) up = String_To_Dbl(n_cur -> value); 

      node_num = io -> tree -> n_root -> num;
      For(j, n_mon)
        {
          if(!strcmp("@root@", mon_list[j])) io -> mcmc -> monitor[node_num] = YES;
        }
      n_cur = XML_Search_Node_Name("appliesto", NO, n_cur -> parent);
      if(n_cur)
        {
          if(n_cur -> attr) tree -> rates -> calib -> proba[node_num] = String_To_Dbl(n_cur -> attr -> value);
          if(!n_cur -> attr && n_cur -> next == NULL)tree -> rates -> calib -> proba[node_num] = 1.;
          if(!n_cur -> attr && n_cur -> next)
            {
              PhyML_Printf("==You either need to provide information about probability with which calibration \n");
              PhyML_Printf("==applies to a node or you need to apply calibartion only to one node. \n");
              PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
              Exit("\n");
            }
          tree -> rates -> calib -> all_applies_to[tree -> rates -> calib -> n_all_applies_to] -> num = node_num; 
          tree -> rates -> calib -> n_all_applies_to++;
          tree -> rates -> calib -> lower = low;
          tree -> rates -> calib -> upper = up;
          PhyML_Printf("\n. .......................................................................");
          PhyML_Printf("\n");
          PhyML_Printf("\n. Node number to which calibration applies to is: [%d]",node_num);
          PhyML_Printf("\n. Clade name: [%s]", "root");
          PhyML_Printf("\n. Lower bound set to: %15f time units.", low);
          PhyML_Printf("\n. Upper bound set to: %15f time units.", up);
          PhyML_Printf("\n. .......................................................................\n");        
        }
    }
  else
    {
      PhyML_Printf("\n==You need to provide calibration information for the root. \n");
      PhyML_Printf("==Root calibration node name in xml file must be 'calibration_root'. \n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
 
  tree -> rates -> calib = last_calib;
  while(tree -> rates -> calib -> prev) tree -> rates -> calib = tree -> rates -> calib -> prev;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //Check for the sum of probabilities for one calibration add up to one
  do
    {
      phydbl p = 0.0;
      for(i = tree -> n_otu; i < 2 * tree -> n_otu - 1; i++) p = p + tree -> rates -> calib -> proba[i];
      if(p != 1.0)
        {
          PhyML_Printf("\n==You need to provide proper probabilities for the calibration. \n");
          PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");
        } 
      if(tree -> rates -> calib -> next) tree -> rates -> calib = tree -> rates -> calib -> next;
      else break;
    }
  while(tree -> rates -> calib);
  while(tree -> rates -> calib -> prev) tree -> rates -> calib = tree -> rates -> calib -> prev;
  ////////////////////////////////////////////////////////////////////////////////////////////////  
  //Set_Current_Calibration(0, tree);
  //TIMES_Set_All_Node_Priors(tree);
  //Slicing_Calibrations(tree);
  //Exit("\n");
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////   
  //clear memory:
  free(clade_name);
  For(i, T_MAX_FILE)
    {
      free(clade[i]);
    }
  free(clade);

  For(i, T_MAX_FILE)
    {
      free(mon_list[i]);
    }
  free(mon_list);

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////   
  //START analysis:
  r_seed = (io -> r_seed < 0)?(time(NULL)):(io -> r_seed);
  srand(r_seed); 
  rand(); 	
  PhyML_Printf("\n. Seed: %d\n", r_seed);
  PhyML_Printf("\n. Pid: %d\n",getpid()); 
  PhyML_Printf("\n. Compressing sequences...\n");
  data = io -> data;
  data[0] -> len = strlen(data[0] -> state); 
  ////////////////////////////////////////////////////////////////////////////
  //memory for compressed sequences:
  cdata         = (calign *)mCalloc(1,sizeof(calign));
  c_seq   = (align **)mCalloc(io -> n_otu,sizeof(align *));
  For(i, io -> n_otu)
    {
      c_seq[i]          = (align *)mCalloc(1,sizeof(align));
      c_seq[i] -> name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      c_seq[i] -> state = (char *)mCalloc(data[0] -> len + 1,sizeof(char));
    }
  cdata -> c_seq = c_seq;
  ////////////////////////////////////////////////////////////////////////////
  cdata = Compact_Data(data, io);
  Free_Seq(io -> data, cdata -> n_otu);
  io -> mod -> io = io;
  Check_Ambiguities(cdata, io -> mod -> io -> datatype, io -> state_len);            
  Make_Model_Complete(mod);                           
  Init_Model(cdata, mod, io);
  if(io -> mod -> use_m4mod) M4_Init_Model(mod -> m4mod, cdata, mod);
  time(&(t_beg));

  tree -> mod         = mod;                                                                    
  tree -> io          = io;
  tree -> data        = cdata;
  tree -> n_pattern   = tree -> data -> crunch_len / tree -> io -> state_len;
  
  Set_Both_Sides(YES, tree);
  Prepare_Tree_For_Lk(tree);
  Set_Current_Calibration(0, tree);
  TIMES_Set_All_Node_Priors(tree);
  
  //set initial value for Hastings ratio for conditional jump:
  tree -> rates -> c_lnL_Hastings_ratio = 0.0;
 
  TIMES_Get_Number_Of_Time_Slices(tree);
  TIMES_Label_Edges_With_Calibration_Intervals(tree);


  tree -> write_br_lens = NO;	

  PhyML_Printf("\n");														
  PhyML_Printf("\n. Input tree with calibration information ('almost' compatible with MCMCtree).\n");
  PhyML_Printf("\n. %s \n", Write_Tree(tree, YES));

  tree -> write_br_lens = YES;


  // Work with log of branch lengths?
  if(tree -> mod -> log_l == YES) Log_Br_Len(tree);

  if(io -> mcmc -> use_data == YES)																
    { 
      // Force the exact likelihood score 
      user_lk_approx = tree -> io -> lk_approx;													
      tree -> io -> lk_approx = EXACT;
		      
      // MLE for branch lengths 																                      
      /* printf("\n. %s",Write_Tree(tree,NO)); */
      /* printf("\n. alpha %f",tree->mod->ras->alpha->v); */
      /* printf("\n.  %f %f %f %f",tree->mod->e_frq->pi->v[0],tree->mod->e_frq->pi->v[1],tree->mod->e_frq->pi->v[2],tree->mod->e_frq->pi->v[3]); */
      /* Lk(NULL,tree); */
      /* printf("\n. %f",tree->c_lnL); */
      /* Exit("\n"); */

      PhyML_Printf("\n");
      Round_Optimize(tree, tree -> data, ROUND_MAX);
 		      
      // Set vector of mean branch lengths for the Normal approximation of the likelihood 
      RATES_Set_Mean_L(tree);
      
      // Estimate the matrix of covariance for the Normal approximation of the likelihood
      PhyML_Printf("\n");
      PhyML_Printf("\n. Computing Hessian...");												    
      tree -> rates -> bl_from_rt = 0;																		
      Free(tree -> rates -> cov_l);																			
      tree -> rates -> cov_l = Hessian_Seo(tree);
															 
      // tree->rates->cov_l = Hessian_Log(tree); 
      For(i, (2 * tree -> n_otu - 3) * (2 * tree -> n_otu - 3)) tree -> rates -> cov_l[i] *= -1.0;
      if(!Iter_Matinv(tree -> rates -> cov_l, 2 * tree -> n_otu - 3, 2 * tree -> n_otu - 3, YES)) Exit("\n");			
      tree -> rates -> covdet = Matrix_Det(tree -> rates -> cov_l, 2 * tree -> n_otu - 3, YES);							
      For(i,(2 * tree -> n_otu - 3) * (2 * tree -> n_otu - 3)) tree -> rates -> invcov[i] = tree -> rates -> cov_l[i];  
      if(!Iter_Matinv(tree -> rates -> invcov, 2 * tree -> n_otu - 3, 2 * tree -> n_otu - 3, YES)) Exit("\n");
      tree -> rates -> grad_l = Gradient(tree);					

      // Pre-calculation of conditional variances to speed up calculations 					
      RATES_Bl_To_Ml(tree);
      RATES_Get_Conditional_Variances(tree);
      RATES_Get_All_Reg_Coeff(tree);
      RATES_Get_Trip_Conditional_Variances(tree);
      RATES_Get_All_Trip_Reg_Coeff(tree);
		      
      Lk(NULL, tree);
      PhyML_Printf("\n");
      PhyML_Printf("\n. p(data|model) [exact ] ~ %.2f",tree -> c_lnL);																			
      tree -> io -> lk_approx = NORMAL;															
      For(i,2 * tree -> n_otu - 3) tree -> rates -> u_cur_l[i] = tree -> rates -> mean_l[i] ;				
      tree -> c_lnL = Lk(NULL,tree);
      PhyML_Printf("\n. p(data|model) [approx] ~ %.2f",tree -> c_lnL);

      tree -> io -> lk_approx = user_lk_approx;	
										
    }

  
  tree -> rates -> model = io -> rates -> model;
  													  
  PhyML_Printf("\n. Selected model '%s' \n", RATES_Get_Model_Name(io -> rates -> model));

  if(tree -> rates -> model == GUINDON) tree -> mod -> gamma_mgf_bl = YES;
														
  tree -> rates -> bl_from_rt = YES;

  if(tree -> io -> cstr_tree) Find_Surviving_Edges_In_Small_Tree(tree, tree -> io -> cstr_tree); 																				 
  time(&t_beg);

  tree -> mcmc = MCMC_Make_MCMC_Struct();
 
  MCMC_Copy_MCMC_Struct(tree -> io -> mcmc, tree -> mcmc, "phytime"); 

  tree -> mod -> m4mod = m4mod;
		
  MCMC_Complete_MCMC(tree -> mcmc, tree);

  tree -> mcmc -> is_burnin = NO;

  //tree -> rates -> nd_t[tree -> n_root -> num] = -1.1;
  //PhyML_Printf("\n");
  //PhyML_Printf("\n. Computing Normalizing Constant(s) for the Node Times Prior Density...\n");
  //tree -> K = Norm_Constant_Prior_Times(tree);
  //Exit("\n");

  MCMC(tree);   

  //PhyML_Printf("\n. Normolizing constant [%f] [%f]\n", tree -> K[0], tree -> K[1]);
                                         															
  MCMC_Close_MCMC(tree -> mcmc);																	
  MCMC_Free_MCMC(tree -> mcmc);														 				
  PhyML_Printf("\n");	
 
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);
  Free_Cseq(cdata);
  Free_Model(mod);
  if(io -> fp_in_align)   fclose(io -> fp_in_align);
  if(io -> fp_in_tree)    fclose(io -> fp_in_tree);
  if(io -> fp_out_lk)     fclose(io -> fp_out_lk);
  if(io -> fp_out_tree)   fclose(io -> fp_out_tree);
  if(io -> fp_out_trees)  fclose(io -> fp_out_trees);
  if(io -> fp_out_stats)  fclose(io -> fp_out_stats);
  fclose(f);
  Free(most_likely_tree);
  Free_Input(io);
  Free_Calib(tree -> rates -> calib);
  time(&t_end);
  Print_Time_Info(t_beg,t_end);	
 	
  /* return 1;    */
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//Function calculates the normalizing constant K of the joint distribution Yule_Order.
//Use the fact that density Yule_Order can be used streight forward only in case of the complete 
//overlap of the calibration intervals for all of the nodes or in case of no overlap.  
phydbl Slicing_Calibrations(t_tree *tree)
{
  int i, j, k, f, n_otu, *indic, *n_slice, *slice_numbers;
  phydbl K, buf, chop_bound, *t_prior_min, *t_prior_max, *t_slice, *t_slice_min, *t_slice_max;


  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  n_otu = tree -> n_otu;

  t_slice        = (phydbl *)mCalloc(2 * (n_otu - 1), sizeof(phydbl)); //the vector of the union of lower and upper bounds, lined up in incresing order.
  t_slice_min    = (phydbl *)mCalloc(2 * n_otu - 3, sizeof(phydbl));   //vector of the lower bounds of the sliced intervals.
  t_slice_max    = (phydbl *)mCalloc(2 * n_otu - 3, sizeof(phydbl));   //vector of the upper bounds of the sliced intervals.
  indic          = (int *)mCalloc((n_otu - 1) * (2 * n_otu - 3), sizeof(int));  //vector of the indicators, columns - node numbers (i + n_otu), rows - the number of the sliced interval.
  slice_numbers  = (int *)mCalloc((n_otu - 1) * (2 * n_otu - 3), sizeof(int )); //vecor of the slice intervals numbers, columns node numbers (i + n_otu), rows - the number of the sliced interval.
  n_slice        = (int *)mCalloc(n_otu - 1, sizeof(int));                      //vector of the numbers of sliced intervals that apply to one node with number (i + n_otu).
  
  i = 0;
  K = 0;
  j = n_otu;
  ////////////////////////////////////////////////////////////////////////////
  //Put prior bounds in one vector t_slice. Excluding tips.
  For(i, n_otu - 1)  
    {
      t_slice[i] = t_prior_min[j];
      j++;
    }

  j = n_otu; 
  for(i = n_otu - 1; i < 2 * n_otu - 3; i++) 
    {
      t_slice[i] = t_prior_max[j];
      j++;
    }
  if(tree -> rates -> nd_t[tree -> n_root -> num] > t_prior_min[tree -> n_root -> num]) chop_bound =  MIN(tree -> rates -> nd_t[tree -> n_root -> num], t_prior_max[tree -> n_root -> num]);
  else chop_bound = t_prior_min[tree -> n_root -> num];
  t_slice[2 * n_otu - 3] = chop_bound; 
  //t_slice[2 * n_otu - 3] = -1.1; 
  //For(j, 2 * n_otu - 2) printf("\n. Slice bound [%f] \n", t_slice[j]);
  ////////////////////////////////////////////////////////////////////////////
  //Get slices in increasing order. Excluding tips.
  do
    {
      f = NO;
      For(j, 2 * n_otu - 3)
        {
          if(t_slice[j] > t_slice[j + 1])
            {
              buf = t_slice[j];
              t_slice[j] = t_slice[j + 1];
              t_slice[j + 1] = buf;
              f = YES;
            }
        }
    }
  while(f);
  //For(j, 2 * n_otu - 2) printf("\n. [1] Slice bound [%f] \n", t_slice[j]);
  for(j = 1; j < 2 * n_otu - 2; j++) t_slice[j] = MAX(chop_bound, t_slice[j]);
  //for(j = 1; j < 2 * n_otu - 2; j++) t_slice[j] = MAX(-1.1, t_slice[j]);
  //For(j, 2 * n_otu - 2) printf("\n. [2] Slice bound [%f] \n", t_slice[j]);
  ////////////////////////////////////////////////////////////////////////////
  //Get the intervals with respect to slices. Total number of t_slice_min(max) - 2 * n_otu - 3. Excluding tips.
  i = 0;
  For(j, 2 * n_otu - 3)
    {
      t_slice_min[j] = t_slice[i];
      t_slice_max[j] = t_slice[i + 1];
      i++;
    }

  //For(j, 2 * n_otu - 3) printf("\n. The interval number [%d] min [%f] max[%f] \n", j, t_slice_min[j], t_slice_max[j]);

  ////////////////////////////////////////////////////////////////////////////
  //Getting indicators for the node number [i + n_otu] to have slice. i = i + n_otu is the node number on the tree and j is the slice number, total 
  //number of intervals is 2 * n_otu - 3. Excluding tips.
  For(i, n_otu - 1)
    { 
      For(j, 2 * n_otu - 3)
        {

          if(Are_Equal(t_prior_min[i + n_otu], t_slice_min[j], 1.E-10) && t_prior_max[i + n_otu] > t_slice_max[j] && t_prior_min[i + n_otu] < t_slice_max[j] && !Are_Equal(t_slice_max[j], t_slice_min[j], 1.E-10)) indic[i * (2 * n_otu - 3) + j] = 1;
          else if(Are_Equal(t_prior_max[i + n_otu], t_slice_max[j], 1.E-10) && t_prior_min[i + n_otu] < t_slice_min[j] && t_prior_max[i + n_otu] > t_slice_min[j] && !Are_Equal(t_slice_max[j], t_slice_min[j], 1.E-10)) indic[i * (2 * n_otu - 3) + j] = 1;
          else if(t_prior_min[i + n_otu] < t_slice_min[j] && t_prior_max[i + n_otu] > t_slice_max[j] && !Are_Equal(t_slice_max[j], t_slice_min[j], 1.E-10)) indic[i * (2 * n_otu - 3) + j] = 1;
          else if(Are_Equal(t_prior_min[i + n_otu], t_slice_min[j], 1.E-10) && Are_Equal(t_prior_max[i + n_otu], t_slice_max[j], 1.E-10)) indic[i * (2 * n_otu - 3) + j] = 1;
        }
    } 

  For(i, n_otu - 2)
    {             
      indic[i * (2 * n_otu - 3)] = 0;
    }
 
  /* 
  printf("\n");
  For(i, n_otu - 1)
    {
      printf(" ['%d]' ", i + n_otu);
      For(j, 2 * n_otu - 3)  
        {          
          printf(". '%d' ", indic[i * (2 * n_otu - 3) + j]);          
        }
      printf("\n");
    }
  */
  ////////////////////////////////////////////////////////////////////////////
  //Get the number of slices that can be applied for each node and the vectors of slice numbers for each node. 
  For(i, n_otu - 1)
    {
      k = 0;
      For(j, 2 * n_otu - 3)
        {
          if(indic[i * (2 * n_otu - 3) + j] == 1) 
            {
              slice_numbers[i * (2 * n_otu - 3) + k] = j; //printf("\n. Node [%d] slice'%d' ", i + n_otu, slice_numbers[i * (2 * n_otu - 3) + j]);
              n_slice[i]++; 
              k++;
            }
        }
      //printf(" Number of slices'%d' \n", n_slice[i]);
    }
  /*
  printf("\n");
  For(i, n_otu - 1)
    {  
      printf(" ['%d]' ", i + n_otu);
       For(j, n_slice[i])  
        {          
          printf(". '%d' ", slice_numbers[i * (2 * n_otu - 3) + j]);          
        }
      printf("\n");
    }
  */

  ////////////////////////////////////////////////////////////////////////////
  //Running through all of the combinations of slices
  int l, tot_num_comb, *cur_slices, *cur_slices_shr, shr_num_slices;
  phydbl P, *t_cur_slice_min, *t_cur_slice_max;
  phydbl k_part;

  tot_num_comb = 1;
  P = 0.0; 
 

  t_cur_slice_min    = (phydbl *)mCalloc(n_otu - 1, sizeof(phydbl));
  t_cur_slice_max    = (phydbl *)mCalloc(n_otu - 1, sizeof(phydbl));
  cur_slices         = (int *)mCalloc(n_otu - 1, sizeof(int)); //the vector of the current slices with repetition.
  cur_slices_shr     = (int *)mCalloc(n_otu - 1, sizeof(int)); //the vector of the current slices without repetition.

  For(i, n_otu - 1) tot_num_comb = tot_num_comb * n_slice[i]; //printf("\n. Total number of combinations of slices [%d] \n", tot_num_comb);
  
  For(k, tot_num_comb)
    {
      shr_num_slices = 0;
      //printf("\n"); 
      For(i, n_otu - 1) //node number i + n_otu
        {
          //printf(" ['%d]' ", i + n_otu);
          l = (k % Number_Of_Comb_Slices(i, n_otu - 1, n_slice)) / Number_Of_Comb_Slices(i+1, n_otu - 1, n_slice); //printf(" Slice number'%d' ", slice_numbers[i * (2 * n_otu - 3) + l]);
          t_cur_slice_min[i] = t_slice_min[slice_numbers[i * (2 * n_otu - 3) + l]]; //printf(" '%f' ", t_cur_slice_min[i]);
          t_cur_slice_max[i] = t_slice_max[slice_numbers[i * (2 * n_otu - 3) + l]]; //printf(" '%f' ", t_cur_slice_max[i]);
          cur_slices[i] = slice_numbers[i * (2 * n_otu - 3) + l];
          //printf("\n"); 
        }
      //printf("\n");
      //For(i, n_otu - 1) printf(" Slice number'%d' ", cur_slices[i]); 
      //printf("\n");

      ///////////////////////////////////////////////////////////////////////////
      //Taking away duplicated slices
      For(i, n_otu - 1)
        {
          for(j = i + 1; j < n_otu - 1; j++)
            {
              if(cur_slices[i] == cur_slices[j]) cur_slices[j] = -1;
            }  
        }
      //For(i, n_otu - 1) printf(" Slice number'%d' \n", cur_slices[i]); 

      ///////////////////////////////////////////////////////////////////////////
      //Getting a vector of all of the slices without duplicates.
      For(i, n_otu -1)
        {
          if(cur_slices[i] >= 0) 
            {
              cur_slices_shr[shr_num_slices] = cur_slices[i];
              shr_num_slices++; 
            }
        }
      //printf("\n"); 
      //For(i, shr_num_slices) printf("\n. Slice number'%d' \n", cur_slices_shr[i]); 
      //printf("\n"); 

      ////////////////////////////////////////////////////////////////////////////
      //Check for the time slices to be set properly
      int result;

      result = TRUE;

      Check_Time_Slices(tree -> n_root, tree -> n_root -> v[1], &result, t_cur_slice_min, t_cur_slice_max, tree);
      Check_Time_Slices(tree -> n_root, tree -> n_root -> v[2], &result, t_cur_slice_min, t_cur_slice_max, tree);
      //printf("\n. '%d' \n", result);

      //For(i, n_otu - 1) printf("\n. Node [%d] min [%f] max [%f] \n", i + n_otu, t_cur_slice_min[i], t_cur_slice_max[i]);

      ////////////////////////////////////////////////////////////////////////////
      //Calculating k_part

      k_part = 1.0;

      if(result != TRUE) k_part = 0.0; 
      else
        {
          int n_1, n_2;

          ////////////////////////////////////////////////////////////////////////////
          //Getting the root node in a given slice
          int *root_nodes;
          int num_elem;
           
          num_elem = 0;
           
          root_nodes = (int *)mCalloc(n_otu - 1, sizeof(int)); 
          
          //printf("\n. Number of slices shrinked [%d] \n", shr_num_slices);
          
          For(i, shr_num_slices)
            {
              //printf("\n. Hello [%d] \n", i);
              //printf("\n. The number of the shrinked interval [%d] min [%f] max [%f] \n", cur_slices_shr[i], t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]]);
              Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[1], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_cur_slice_min, t_cur_slice_max, tree);
              Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[2], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_cur_slice_min, t_cur_slice_max, tree); 
            }
          //printf("\n. Number of elements in a vector of the nodes [%d] \n", num_elem);
          //For(j, num_elem) printf("\n. Root node number [%d] \n", root_nodes[j] + n_otu);
          For(j, num_elem) 
            {
              //printf("\n. Root node number [%d] \n", tree -> a_nodes[root_nodes[j] + n_otu] -> num); 
              
              n_1 = 0;
              n_2 = 0;
              Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[1], &n_1, t_cur_slice_min, t_cur_slice_max, tree);
              Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[2], &n_2, t_cur_slice_min, t_cur_slice_max, tree);
              //printf("\n. n_1 [%d] n_2 [%d]\n", n_1, n_2);  
              k_part = k_part * Factorial(n_1 + n_2) / ((phydbl)Factorial(n_1 + n_2 + 1) * Factorial(n_1) * Factorial(n_2));
            }
          //printf("\n. k_part [%f] \n", k_part);       
          ////////////////////////////////////////////////////////////////////////////
          //Calculating PRODUCT over all of the time slices k_part * (exp(-lmbd*l) - exp(-lmbd*u))/(exp(-lmbd*l) - exp(-lmbd*u)) 
          phydbl num, denom, lmbd;
          
          lmbd = tree -> rates -> birth_rate;
          num = 1;
          denom = 1;
          //lmbd = 6.001;
          /* For(j, n_otu - 1) num = num * (EXP(-lmbd * t_cur_slice_min[j]) - EXP(-lmbd * t_cur_slice_max[j]));  */
          /* for(j = n_otu; j < 2 * n_otu - 1; j++) denom = denom * (EXP(-lmbd * t_prior_min[j]) - EXP(-lmbd * t_prior_max[j]));  */
          For(j, n_otu - 2) num = num * (EXP(lmbd * t_cur_slice_max[j]) - EXP(lmbd * t_cur_slice_min[j])); 
          for(j = n_otu; j < 2 * n_otu - 2; j++) denom = denom * (EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j])); 
          k_part = (k_part * num) / denom; 
          //printf("\n. [2] k_part of the tree for one combination of slices [%f] \n", k_part); 
        } 
      P = P + k_part;     
    }
  //printf("\n. [P] of the tree for one combination of slices [%f] \n", P);
  K = 1 / P;
  //printf("\n. [K] of the tree for one combination of slices [%f] \n", K);
  ////////////////////////////////////////////////////////////////////////////
  free(t_cur_slice_min);
  free(t_cur_slice_max);
  free(cur_slices);
  free(cur_slices_shr); 
  free(t_slice);        
  free(t_slice_min);    
  free(t_slice_max);    
  free(indic);          
  free(slice_numbers);  
  free(n_slice);
  //free(tree -> K);
        
  return(K); 
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int Number_Of_Comb_Slices(int m, int num_elem, int *n_slice)
{
  int i, num_comb;

  i = 0;
  num_comb = 1;       

  for(i = m; i < num_elem; i++) num_comb = num_comb * n_slice[i]; 
  
  return(num_comb);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Check the combination of the time slices to be set correctly.
void Check_Time_Slices(t_node *a, t_node *d, int *result, phydbl *t_cur_slice_min, phydbl *t_cur_slice_max, t_tree *tree)
{
  int n_otu;

  n_otu = tree -> n_otu;
 

  d -> anc = a;
  if(d -> tax) return;
  else
    { 
      if(t_cur_slice_max[d -> num - n_otu] < t_cur_slice_max[a -> num - n_otu])
        {
          *result = FALSE; 
        }

      int i;
      For(i,3) 
	if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root))
             Check_Time_Slices(d, d -> v[i], result, t_cur_slice_min, t_cur_slice_max, tree);
    }
}

//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//Getting the number of nodes on both sides from the node d_start, that are in the slice of that node.
void Number_Of_Nodes_In_Slice(t_node *d_start, t_node *d, int *n, phydbl *t_cur_slice_min, phydbl *t_cur_slice_max, t_tree *tree)
{
  int n_otu;
 
  n_otu = tree -> n_otu;


  if(d -> tax) return;
  else
    { 
      if(Are_Equal(t_cur_slice_max[d_start -> num - n_otu], t_cur_slice_max[d -> num - n_otu], 1.E-10) && Are_Equal(t_cur_slice_min[d_start -> num - n_otu], t_cur_slice_min[d -> num - n_otu], 1.E-10))
        {
          (*n)++; 
          int i;
          For(i,3) 
            if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root))
              Number_Of_Nodes_In_Slice(d_start, d -> v[i], n, t_cur_slice_min, t_cur_slice_max, tree);
        }      
    }
}


//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//Returnig the root node in the given slice.
void Search_Root_Node_In_Slice(t_node *d_start, t_node *d, int *root_nodes, int *num_elem, phydbl t_slice_min, phydbl t_slice_max, phydbl *t_cur_slice_min, phydbl *t_cur_slice_max, t_tree *tree)
{
  int j, n_otu, f;
  j = 0;
  f = FALSE;
  n_otu = tree -> n_otu; //printf("\n. Node number 1  [%d] \n", d_start -> num);printf("\n. Node number 1  [%d] \n", d -> num);

  if(Are_Equal(t_cur_slice_max[d_start -> num - n_otu], t_slice_max, 1.E-10) && Are_Equal(t_cur_slice_min[d_start -> num - n_otu], t_slice_min, 1.E-10))
    {
      For(j, *num_elem) if(d_start -> num - n_otu == (root_nodes[j])) f = TRUE;
      if(f != TRUE)
        {
          (root_nodes[(*num_elem)]) = d_start -> num - n_otu; //printf("\n. Node number 2_2  [%d] \n", root_nodes[(*num_elem)] + n_otu);
          (*num_elem)++;  //printf("\n. Number of elements 2_2  [%d] \n", *num_elem);
          return;
        }
      
    }
  else
    {
      d -> anc = d_start;
      if(d -> tax) return;
      else
        { 
          if(Are_Equal(t_cur_slice_max[d -> num - n_otu], t_slice_max, 1.E-10) && Are_Equal(t_cur_slice_min[d -> num - n_otu], t_slice_min, 1.E-10))
            {
              (root_nodes[*num_elem]) = d -> num - n_otu; //printf("\n. Node number 3_2  [%d] \n",  root_nodes[*num_elem] + n_otu);
              (*num_elem)++; //printf("\n. Number of elements 3_2  [%d] \n", *num_elem);  
              return; 
            }
          
          int i;
          For(i,3) 
            if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root))
              Search_Root_Node_In_Slice(d, d -> v[i], root_nodes, num_elem, t_slice_min, t_slice_max, t_cur_slice_min, t_cur_slice_max, tree);
        }
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int Factorial(int base)
{
  if(base == 0) return(1); 
  if(base == 1) return(1);
  return(base * Factorial(base-1)); 
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Calculate the vector of the norm.constants for prior on node times. 
//The length of the vector is the total number of combinations of calibrations.
phydbl *Norm_Constant_Prior_Times(t_tree *tree)
{

  phydbl *t_prior_min, *t_prior_max, *K; 
  short int *t_has_prior;
  int i, j, k, tot_num_comb;
  t_cal *calib;
 



  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;

  tot_num_comb = Number_Of_Comb(calib);
  
  //PhyML_Printf("\n. The total number of calibration combinations: [%d]\n", tot_num_comb);

  K = (phydbl *)mCalloc(tot_num_comb, sizeof(phydbl));   

  For(i, tot_num_comb)
    {
      for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
        {
          t_prior_min[j] = -BIG;
          t_prior_max[j] = BIG;
          t_has_prior[j] = NO; 
        }
      do
        {
          k = (i % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);
          t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower);
          t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper);
          t_has_prior[calib -> all_applies_to[k] -> num] = YES;
          if(calib -> next) calib = calib -> next;
          else break;
        }
      while(calib);
      TIMES_Set_All_Node_Priors(tree);

      K[i] = Slicing_Calibrations(tree);     
      //PhyML_Printf("\n. Number [%d] normolizing constant [%f] \n", i+1, K[i]);
      while(calib -> prev) calib = calib -> prev;
    }
  return(K);
}
