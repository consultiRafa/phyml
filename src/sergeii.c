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
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Check_Node_Time(t_tree *tree)
{
  phydbl *t_prior_min, *t_prior_max, *nd_t; 
  int i, result;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  nd_t = tree -> rates -> nd_t;
  For(i, 2 * tree -> n_otu - 1)
    {
      if(nd_t[i] >= t_prior_min[i] && nd_t[i] <= t_prior_max[i])
        {
          result = TRUE;
          //printf("\n. node number:'%d' time cur:'%f' time min:'%f' time max:'%f'", i, tree -> rates -> nd_t[i], tree -> rates -> t_prior_min[i], tree -> rates -> t_prior_max[i]);
          //PhyML_Printf("\n\n==You have a mistake in calibration information.\n");
          //PhyML_Printf("\n==Err in file %s at line %d\n",__FILE__,__LINE__);
          //Exit("\n");
        }
      else
        { 
          result = FALSE;
          //printf("\n. '%d' \n", i);
          break;
        }
    }
  return(result);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

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

void Set_Current_Calibration(int row, t_tree *tree)
{

  t_cal *calib;
  phydbl *t_prior_min, *t_prior_max; 
  short int *t_has_prior;
  int k;

  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  k = -1;
  do
    {
      k = (row % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);      
      t_prior_min[calib -> all_applies_to[k] -> num] = calib -> lower;
      t_prior_max[calib -> all_applies_to[k] -> num] = calib -> upper;
      t_has_prior[calib -> all_applies_to[k] -> num] = YES; 
      /*PhyML_Printf("\n. .......................................................................");
      PhyML_Printf("\n");
      PhyML_Printf("\n. Node number to which calibration applies to is: [%d]", calib -> all_applies_to[k] -> num);
      PhyML_Printf("\n. Lower bound set to: %15f time units.", t_prior_min[calib -> all_applies_to[k] -> num]);
      PhyML_Printf("\n. Upper bound set to: %15f time units.", t_prior_max[calib -> all_applies_to[k] -> num]);
      PhyML_Printf("\n. .......................................................................");*/
      if(calib->next) calib = calib->next;
      else break;    
    }
  while(calib);
  //while(calib -> prev) calib = calib -> prev;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl TIMES_Calib_Cond_Prob(t_tree *tree)
{

  phydbl times_partial_proba, times_tot_proba, *t_prior_min, *t_prior_max; 
  short int *t_has_prior;
  int i, j, k, tot_num_comb;
  t_cal *calib;
 

  times_partial_proba = 1;
  //times_tot_proba = 0;
  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  tot_num_comb = Number_Of_Comb(calib);

  For(i, tot_num_comb)
    {
      for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
        {
          t_prior_min[j] = -BIG;
          t_prior_max[j] = BIG;
          t_has_prior[j] = NO; 
        }
      times_partial_proba = 1; 
      do
        {
          k = (i % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);
          t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower);
          t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper);
          t_has_prior[calib -> all_applies_to[k] -> num] = YES;
          if((t_prior_min[calib -> all_applies_to[k] -> num] - t_prior_max[calib -> all_applies_to[k] -> num]) > 0) times_partial_proba = 0; 
          else times_partial_proba *= calib -> proba[calib -> all_applies_to[k] -> num]; 
          if(calib -> next) calib = calib -> next;
          else break;
        }
      while(calib);
      TIMES_Set_All_Node_Priors(tree); 

      if(Check_Node_Time(tree) != TRUE) times_partial_proba = 0;      
      /* if(Check_Calibration_Consistency(tree) != TRUE) times_partial_proba = .0; */
      
      times_tot_proba += (times_partial_proba * EXP(TIMES_Lk_Yule_Order(tree)));
      
      while(calib -> prev) calib = calib -> prev;
    }

  Set_Current_Calibration(0, tree);   
  TIMES_Set_All_Node_Priors(tree);     

  return(LOG(times_tot_proba));
}


void PhyTime_XML(char *xml_file)
{

  FILE *f;
  char **clade, *clade_name, **mon_list;
  phydbl low, up; //times_tot_proba;
  int i, j, n_taxa, clade_size, node_num, n_mon, tot_num_cal;
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
 

  i = 0;
  j = 0;
  tot_num_cal = 0;
  last_calib       = NULL;
  mod              = NULL;
  most_likely_tree = NULL;
   

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
  Update_Ancestors(io -> tree -> n_root, io -> tree -> n_root -> v[0], io -> tree);
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

  n_r = n_r -> child; 

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
          tot_num_cal++;
	  if (tree -> rates -> calib == NULL) tree -> rates -> calib = Make_Calib(tree -> n_otu);
          if(last_calib)
            {
              last_calib -> next = tree -> rates -> calib;
              tree -> rates -> calib -> prev = last_calib;
            }
          last_calib = tree -> rates -> calib;

	  low = BIG;
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
                  //PhyML_Printf("\n. '%f'\n", tree -> rates -> calib -> proba[node_num]);
                  tree -> rates -> calib -> all_applies_to[tree -> rates -> calib -> n_all_applies_to] -> num = node_num; 
                  //PhyML_Printf("\n. '%d'\n", tree -> rates -> calib -> all_applies_to[tree -> rates -> calib -> n_all_applies_to] -> num);
                  tree -> rates -> calib -> n_all_applies_to++;
                  tree -> rates -> calib -> lower = low;
                  //PhyML_Printf("\n. '%f'\n", tree -> rates -> calib -> lower);
                  tree -> rates -> calib -> upper = up;
                  //PhyML_Printf("\n. '%f'\n", tree -> rates -> calib -> upper);
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
                  //printf(" '%f' '%f' \n", low, up); 
                  //printf(" '%d' \n", node_num);   
                  //tree -> rates -> t_prior_min[node_num] = low;
                  //tree -> rates -> t_prior_max[node_num] = up;
                  //tree -> rates -> t_has_prior[node_num] = YES;
                  //printf(" '%f' '%f' \n", tree -> rates -> t_prior_min[node_num], tree -> rates -> t_prior_max[node_num]);              
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
 
  tree -> rates -> calib = last_calib;
  while(tree -> rates -> calib -> prev) tree -> rates -> calib = tree -> rates -> calib -> prev;

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////   
  /*Set_Current_Calibration(1, tree);
  TIMES_Set_All_Node_Priors(tree);
  MCMC_Randomize_Node_Times(tree); 
  tree -> rates -> t_prior_min[109] = -10;
  tree -> rates -> t_prior_max[109] = 0;
  tree -> rates -> nd_t[109] = 1;
  for(j = tree -> n_otu; j < 2*tree->n_otu-1; j++) printf("\n. node number:'%d' time:'%f'", j, tree -> rates -> nd_t[j]);
  printf("\n. '%d' \n", Check_Node_Time(tree));
  Exit("\n");
  for(j = tree -> n_otu; j < 2*tree->n_otu-1; j++) printf("\n. node number:'%d' time min:'%f' time max:'%f'", j, tree -> rates -> t_prior_min[j], tree -> rates -> t_prior_max[j]); 
  printf("\n\n");
  Set_Current_Calibration(2, tree);
  TIMES_Set_All_Node_Priors(tree); 
  for(j = tree -> n_otu; j < 2*tree->n_otu-1; j++) printf("\n. node number:'%d' time min:'%f' time max:'%f'", j, tree -> rates -> t_prior_min[j], tree -> rates -> t_prior_max[j]); 
  Exit("\n");*/

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

  MCMC(tree);                                            															
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
