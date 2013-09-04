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
  node_num = -1;

  //file can/cannot be open:
  if ((f =(FILE *)fopen(xml_file, "r")) == NULL)
    {
      PhyML_Printf("\n== File '%s' can not be opened...\n",xml_file);
      Exit("\n");
    }

  n_r = XML_Load_File(f);

  //XML_Search_Node_Attribute_Value_Clade("id", "clade1", NO, n_r -> child);
  //Exit("\n");

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
  /* tree        = (t_tree *)mCalloc(1,sizeof(t_tree)); */
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
	  PhyML_Printf("\n== No input tree was not found. \n");
	  PhyML_Printf("\n== Please specify either a tree file name or enter the whole tree directly in the XML file (in Newick format). \n");
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
   
  //TO DO: check that the tree has a root
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
  /* clade = (char **)mCalloc(tree -> n_otu, sizeof(char *)); */
  /* For(i, tree -> n_otu) */
  /*   { */
  /*     clade[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
  /*   } */

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
      PhyML_Printf("\n== There is no calibration information provided. \n");
      PhyML_Printf("\n== Please check your data. \n");
      Exit("\n");
    }
  else
    {
      if(XML_Search_Node_Name("upper", NO, n_cur -> child) == NULL && XML_Search_Node_Name("lower", NO, n_cur -> child) == NULL)
	{
	  PhyML_Printf("\n== There is no calibration information provided. \n");
	  PhyML_Printf("\n== Please check your data. \n");
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
		PhyML_Printf("\n== Sequence type (nt / aa) must be provided. \n");
		PhyML_Printf("\n== Please, include this information as an attribute of component <%s>. \n", n_r -> name);
		Exit("\n");
	  }	
          char *s;
          s = To_Upper_String(n_r -> attr -> value);
	  /* if(!strcmp(To_Upper_String(n_r -> attr -> value), "NT")) */
          if(!strcmp(s, "NT"))
	    {
	      io -> datatype = 0;
	      io -> mod -> ns = 4;
	    }
 	  /* if(!strcmp(To_Upper_String(n_r -> attr -> value), "AA")) */
          if(!strcmp(s, "AA"))
	    {
	      io -> datatype = 1;
	      io -> mod -> ns = 20;
	    }
          free(s);
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
              char *s;	
	      i = 0;
              s = To_Upper_String(n_r -> child -> value);
	      do
		{
                  strcpy(io -> in_align_file, "sergeii"); 
		  strcpy(io -> data[i] -> name, n_r -> child -> attr -> value);
		  /* strcpy(io -> data[i] -> state, To_Upper_String(n_r -> child -> value)); */
                  strcpy(io -> data[i] -> state, s);
		  i++;
		  if(n_r -> child -> next) n_r -> child = n_r -> child -> next;
		  else break;
		}
	      while(n_r -> child);
	      n_taxa = i;
              free(s);
	    }

	  //checking if a sequences of the same lengths:

	  i = 1;
	  For(i, n_taxa) if(strlen(io -> data[0] -> state) != strlen(io -> data[i] -> state))
	    {
	      PhyML_Printf("\n== All the sequences should be of the same length. Please check your data...\n");
	      Exit("\n");
	      break;
	    }  

	  //checking sequence names:
	  i = 0;
	  For(i, n_taxa) Check_Sequence_Name(io -> data[i] -> name);

	  //check if a number of tips is equal to a number of taxa:
	  if(n_taxa != io -> n_otu)
	    {
	      PhyML_Printf("\n== The number of taxa is not the same as a number of tips. Check your data...\n");
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
                  //case of internal node:
                  strcpy(clade_name, n_r -> child -> attr -> value);//reached clade names n_r -> child -> attr -> value
                  if(!strcmp("root", clade_name))
                    {
                      node_num = io -> tree -> n_root -> num;
                      //printf("\n. Node number [%d] \n", node_num);
                    }
                  else if(strcmp("NO_CLADE", clade_name) && strcmp("root", clade_name))
                    {
                      xml_node *n_clade, *nd2;
                      nd2 = n_r -> parent;
                      /* n_clade = XML_Search_Node_Attribute_Value_Clade("id", clade_name, NO, n_r -> parent -> child); */
                      n_clade = XML_Search_Node_Generic("clade", "id", clade_name, YES, nd2);
                      if(n_clade) //found clade with a given name
                        {
                          i = 0;
                          xml_node *nd;
                          nd = n_clade -> child;
                          /* do */
                          /*   { */
                          /*     strcpy(clade[i], nd -> attr -> value);  */
                          /*     i++; */
                          /*     nd = nd -> next; */
                          /*     if(!nd) break; */
                          /*   } */
                          /* while(n_clade -> child); */
                          /* clade_size = i; */
                          clade = XML_Reading_Clade(nd, tree);
                          clade_size = XML_Number_Of_Taxa_In_Clade(nd);
                          node_num = Find_Clade(clade, clade_size, io -> tree);
                          /* printf("\n. Clade size [%d] \n", clade_size); */
                          /* printf("\n. Node number [%d] \n", node_num); */
                         }
                      else
                        {
                          PhyML_Printf("\n== The calibration information for clade [%s] was not found.", clade_name);
                          PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
                          Exit("\n");
                        }
                    }
                  if(strcmp("NO_CLADE", clade_name))
                    {
                      For(j, n_mon)
                        {
                          if(!strcmp(clade_name, mon_list[j])) io -> mcmc -> monitor[node_num] = YES;
                        }
                    }
                  /* For(i, clade_size) PhyML_Printf("\n. Clade name [%s] Taxon name: [%s]", clade_name, clade[i]); */
                  if(strcmp("NO_CLADE", clade_name)) 
                    {
                      tree -> rates -> calib -> proba[node_num] = String_To_Dbl(n_r -> child -> attr -> next -> value);
                      
                      if(!n_r -> child -> attr -> next && n_r -> child -> next ==  NULL) tree -> rates -> calib -> proba[node_num] = 1.;
                      if(!n_r -> child -> attr -> next && n_r -> child -> next)
                        {
                          PhyML_Printf("\n== You either need to provide information about the probability with which calibration ");
                          PhyML_Printf("\n== applies to a node or you need to apply calibration only to one node.");
                          PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
                          Exit("\n");
                        }
                      
                      tree -> rates -> calib -> all_applies_to[tree -> rates -> calib -> n_all_applies_to] -> num = node_num;
                    }
                  else tree -> rates -> calib -> proba[2 * tree -> n_otu - 1] = String_To_Dbl(n_r -> child -> attr -> next -> value);
                  tree -> rates -> calib -> n_all_applies_to++;
                  tree -> rates -> calib -> lower = low;
                  tree -> rates -> calib -> upper = up;
                  /* printf("\n. Porbability [%f] \n", String_To_Dbl(n_r -> child -> attr -> next -> value)); */
                  
                  /////////////////////////////////////////////////////////////////////////////////////////////////////
                  PhyML_Printf("\n. .......................................................................");
                  PhyML_Printf("\n");
                  if(strcmp(clade_name, "NO_CLADE")) PhyML_Printf("\n. Clade name: [%s]", clade_name);
                  else 
                    {
                      PhyML_Printf("\n. Calibration with:");
                      PhyML_Printf("\n. Lower bound set to: %15f time units.", low);
                      PhyML_Printf("\n. Upper bound set to: %15f time units.", up);
                      PhyML_Printf("\n. does not apply to any node with probability [%f]", String_To_Dbl(n_r -> child -> attr -> next -> value));
                      PhyML_Printf("\n. .......................................................................");
                    }
       
                  if(strcmp(clade_name, "root") && strcmp(clade_name, "NO_CLADE"))
                    {
                      For(i, clade_size) PhyML_Printf("\n. Taxon name: [%s]", clade[i]);
                    }
                  if(strcmp(clade_name, "NO_CLADE")) 
                    {
                      PhyML_Printf("\n. Node number to which calibration applies to is [%d] with probability [%f]", node_num, String_To_Dbl(n_r -> child -> attr -> next -> value));
                      PhyML_Printf("\n. Lower bound set to: %15f time units.", low);
                      PhyML_Printf("\n. Upper bound set to: %15f time units.", up);
                      PhyML_Printf("\n. .......................................................................");
                    }
                  /////////////////////////////////////////////////////////////////////////////////////////////////////
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
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //Check for the sum of probabilities for one calibration add up to one
  do
    {
      phydbl p = 0.0;
      for(i = tree -> n_otu; i < 2 * tree -> n_otu; i++) 
        {
          p = p + tree -> rates -> calib -> proba[i];
          /* PhyML_Printf("\n. # applies to %d \n", tree -> rates -> calib -> n_all_applies_to); */
          /* PhyML_Printf("\n. %f \n", tree -> rates -> calib -> proba[i]); */
        }
      if(!Are_Equal(p, 1.0, 1.E-10))
        {
          PhyML_Printf("\n. .......................................................................");
          PhyML_Printf("\n. WARNING! The sum of the probabilities for the calibration with:");
          PhyML_Printf("\n. Lower bound set to: %15f time units.", tree -> rates -> calib -> lower);
          PhyML_Printf("\n. Upper bound set to: %15f time units.", tree -> rates -> calib -> upper);
          PhyML_Printf("\n. is not equal to 1.");
          PhyML_Printf("\n. The probability that this calibration does not apply to any node is set to [%f].", 1.0 - p);
          PhyML_Printf("\n. .......................................................................");
          tree -> rates -> calib -> proba[2 * tree -> n_otu - 1] = 1.0 - p;
          tree -> rates -> calib -> n_all_applies_to++;
          /* PhyML_Printf("\n==You need to provide proper probabilities for the calibration. \n"); */
          /* PhyML_Printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__); */
          /* Exit("\n"); */
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
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////   
  //clear memory:
  free(clade_name);
  For(i, tree -> n_otu)
    {
      free(clade[i]);
    }
  free(clade);
  For(i, T_MAX_FILE)
    {
      free(mon_list[i]);
    }
  free(mon_list);


  //Exit("\n");
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

  RATES_Fill_Lca_Table(tree);

  tree -> mod         = mod;                                                                    
  tree -> io          = io;
  tree -> data        = cdata;
  tree -> n_pattern   = tree -> data -> crunch_len / tree -> io -> state_len;
  
  Set_Both_Sides(YES, tree);
  Prepare_Tree_For_Lk(tree);


  //calculate the probabilities of each combination of calibrations:
  TIMES_Calib_Partial_Proba(tree);
  int cal_numb = 0;
  do
    {
      if(!Are_Equal(tree -> rates -> times_partial_proba[cal_numb], 0.0, 1.E-10)) break;
      else cal_numb += 1;
    }
  while(1);  
  /* printf("\n. Calib number [%d] \n", cal_numb); */
  Set_Current_Calibration(cal_numb, tree);
  int tot_num_comb;
  tot_num_comb = Number_Of_Comb(tree -> rates -> calib);														
  PhyML_Printf("\n. The number of calibration combinations that is going to be considered is %d.\n", tot_num_comb);
  TIMES_Set_All_Node_Priors(tree);

  //set initial value for Hastings ratio for conditional jump:
  /* tree -> rates -> c_lnL_Hastings_ratio = 0.0; */
 
  TIMES_Get_Number_Of_Time_Slices(tree);
  TIMES_Label_Edges_With_Calibration_Intervals(tree);

  tree -> write_br_lens = NO;	
									
  PhyML_Printf("\n. Input tree with calibration information ('almost' compatible with MCMCtree).\n");

  char *t;
  t =  Write_Tree(tree, YES);
  PhyML_Printf("\n. %s \n", t);
  free(t);

  tree -> write_br_lens = YES;


  // Work with log of branch lengths?
  if(tree -> mod -> log_l == YES) Log_Br_Len(tree);

  if(io -> mcmc -> use_data == YES)																
    { 
      // Force the exact likelihood score 
      user_lk_approx = tree -> io -> lk_approx;													
      tree -> io -> lk_approx = EXACT;
		      
      /* MLE for branch lengths 																                       */
      /* printf("\n. %s",Write_Tree(tree,NO)); */
      /* printf("\n. alpha %f",tree->mod->ras->alpha->v); */
      /* printf("\n.  %f %f %f %f",tree->mod->e_frq->pi->v[0],tree->mod->e_frq->pi->v[1],tree->mod->e_frq->pi->v[2],tree->mod->e_frq->pi->v[3]); */
      /* Lk(NULL,tree); */
      /* printf("\n. %f",tree->c_lnL); */
      /* Exit("\n"); */

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

  //PhyML_Printf("\n");
  //PhyML_Printf("\n. Computing Normalizing Constant(s) for the Node Times Prior Density...\n");
  /* tree -> K = Norm_Constant_Prior_Times(tree); */
  /* Exit("\n"); */


  MCMC(tree);   
                                         															
  MCMC_Close_MCMC(tree -> mcmc);																	
  MCMC_Free_MCMC(tree -> mcmc);														 				
  PhyML_Printf("\n");	
 
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);
  Free_Cseq(cdata);
  free(cdata);
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


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Calculate the prior probability for node times taking into account the 
//probailitis with which each calibration applies to the particular node.
phydbl TIMES_Calib_Cond_Prob(t_tree *tree)
{

  phydbl times_lk, *Yule_val, *times_partial_proba, times_tot_proba, *t_prior_min, *t_prior_max, c, constant, ln_t; 
  short int *t_has_prior;
  int i, j, k, tot_num_comb;
  t_cal *calib;
 

  times_tot_proba = 0.0;
  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  times_partial_proba = tree -> rates -> times_partial_proba;
  /* constant = tree -> K; */

  tot_num_comb = Number_Of_Comb(calib);

  
  Yule_val = (phydbl *)mCalloc(tot_num_comb, sizeof(phydbl));    
  /* times_partial_proba = (phydbl *)mCalloc(tot_num_comb, sizeof(phydbl));    */

  For(i, tot_num_comb)
    {
      /* printf(" \n\n. COMBINATION OF CALIBRATIONS [%d] \n\n", i + 1); */
      for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
        {
          t_prior_min[j] = -BIG;
          t_prior_max[j] = BIG;
          t_has_prior[j] = NO; 
        }
      do
        {
          k = (i % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);
          if(calib -> all_applies_to[k] -> num)
            {
              t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower); //printf("min %f",  t_prior_min[calib -> all_applies_to[k] -> num]);
              t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper); //printf("max %f",  t_prior_max[calib -> all_applies_to[k] -> num]);
              t_has_prior[calib -> all_applies_to[k] -> num] = YES;
              /* if((t_prior_min[calib -> all_applies_to[k] -> num] > t_prior_max[calib -> all_applies_to[k] -> num])) times_partial_proba[i] = 0.0;  */
              /* else times_partial_proba[i] *= calib -> proba[calib -> all_applies_to[k] -> num];  */
            }
          /* else */
          /*   { */
          /*     if(calib -> next) calib = calib -> next; */
          /*     else break; */
          /*   } */
          if(calib -> next) calib = calib -> next;
          else break;
        }
      while(calib);

      int result;
      result = TRUE;
      TIMES_Set_All_Node_Priors_S(&result, tree);
      /* TIMES_Set_All_Node_Priors(tree); */
      /* Check_Node_Time(tree -> n_root, tree -> n_root -> v[2], &result, tree); */
      /* Check_Node_Time(tree -> n_root, tree -> n_root -> v[1], &result, tree); */
      /* printf("\n\n"); */
      /* for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) printf("\n. [1] Node [%d] min [%f] max [%f] node time [%f]\n", j, tree -> rates -> t_prior_min[j], tree -> rates -> t_prior_max[j], tree -> rates -> nd_t[j]); */
      /* printf("\n. result = [%d] \n", result); */
      /* printf("\n\n"); */

      /* tree -> rates -> birth_rate = 4.0; */
      /* times_lk = TIMES_Lk_Yule_Order(tree); */
      times_lk =  TIMES_Lk_Yule_Order_Root_Cond(tree);
      /* printf("\n. times_lk = [%f] \n", times_lk); */
      /* printf("\n. result = [%d] \n", result); */
      /* if(result != FALSE) times_lk = TIMES_Lk_Yule_Order(tree); */
      /* else times_lk = 1.0; */

      constant = 1.0; 
      if (tot_num_comb > 1) if(times_lk > -INFINITY && result != FALSE) constant = Slicing_Calibrations(tree);
      /* else */
      /*   { */
      /*     times_lk = 0.0; */
      /*     times_partial_proba[i] = 0.0; */
      /*   } */
      /* printf("\n. p[%i] = %f \n", i + 1, times_partial_proba[i]); */
      /* printf("\n. K = [%f] \n", constant); */
      /* K = Norm_Constant_Prior_Times(tree); */
      /* Yule_val[i] = K[i] * TIMES_Lk_Yule_Order(tree); */
 
      /* For(j, 2 * tree -> n_otu - 1) printf("\n. [2] Node [%d] time [%f]\n", j, tree -> rates -> nd_t[j]); */
 
      Yule_val[i] = LOG(constant) + times_lk;

      /* printf("\n. Yule = %f \n", Yule_val[i]); */
 
      while(calib -> prev) calib = calib -> prev;
      /* Exit("\n"); */
    }
 
  /* min_value = 0.0; */
  /* For(i, tot_num_comb) if(Yule_val[i] < min_value && Yule_val[i] > -INFINITY) min_value = Yule_val[i]; */
  /* c = -600. - min_value;   */

  /* Exit("\n"); */
  c = .0;
  times_tot_proba = 0.0;
  For(i, tot_num_comb)
    {
      times_tot_proba += times_partial_proba[i] * EXP(Yule_val[i] + c);
      /* printf("\n. Proba = [%f] \n", times_partial_proba[i]); */
    } 

  For(i, 2 * tree -> n_otu - 1) t_has_prior[i] = NO;

  ln_t = -c + LOG(times_tot_proba);
  /* printf("\n. Prior for node times = [%f] \n", ln_t); */
  /* Set_Current_Calibration(1, tree); */
  /* printf("\n\n"); */
  free(Yule_val);  
  /* free(times_partial_proba); */
  /* Exit("\n"); */
  return(ln_t);
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//Function calculates the normalizing constant K of the joint distribution Yule_Order.
//Use the fact that density Yule_Order can be used streight forward only in case of the complete 
//overlap of the calibration intervals for all of the nodes or in case of no overlap.  
phydbl Slicing_Calibrations(t_tree *tree)
{
  int i, j, k, f, n_otu, *indic, *n_slice, *slice_numbers;
  phydbl K, buf, chop_bound, *t_prior_min, *t_prior_max, *t_slice, *t_slice_min, *t_slice_max, *g_i_node;


  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  n_otu = tree -> n_otu;

  t_slice        = (phydbl *)mCalloc(2 * (n_otu - 1), sizeof(phydbl)); //the vector of the union of lower and upper bounds, lined up in incresing order.
  t_slice_min    = (phydbl *)mCalloc(2 * n_otu - 3, sizeof(phydbl));   //vector of the lower bounds of the sliced intervals.
  t_slice_max    = (phydbl *)mCalloc(2 * n_otu - 3, sizeof(phydbl));   //vector of the upper bounds of the sliced intervals.
  indic          = (int *)mCalloc((n_otu - 1) * (2 * n_otu - 3), sizeof(int));  //vector of the indicators, columns - node numbers (i + n_otu), rows - the number of the sliced interval.
  slice_numbers  = (int *)mCalloc((n_otu - 1) * (2 * n_otu - 3), sizeof(int )); //vecor of the slice intervals numbers, columns node numbers (i + n_otu), rows - the number of the sliced interval.
  n_slice        = (int *)mCalloc(n_otu - 1, sizeof(int));                      //vector of the numbers of sliced intervals that apply to one node with number (i + n_otu).
  g_i_node       = (phydbl *)mCalloc((n_otu - 1) * (2 * n_otu - 3), sizeof(phydbl));  //vector of the values of the function g evaluated for each node.
  
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
  /* For(i, 2 * n_otu - 2) if(Are_Equal(t_prior_min[tree -> n_root -> num], t_slice[i], 1.E-10)) */
  /* printf("\n. Chop bound [%f] \n", chop_bound);  */
  //t_slice[2 * n_otu - 3] = -1.1; 
  /* For(j, 2 * n_otu - 2) printf("\n. Slice bound [%f] \n", t_slice[j]); */
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
  /* For(j, 2 * n_otu - 2) printf("\n. [1] Slice bound [%f] \n", t_slice[j]); */
  for(j = 1; j < 2 * n_otu - 2; j++) t_slice[j] = MAX(chop_bound, t_slice[j]);
  //for(j = 1; j < 2 * n_otu - 2; j++) t_slice[j] = MAX(-1.1, t_slice[j]);
  /* For(j, 2 * n_otu - 2) printf("\n. [2] Slice bound [%f] \n", t_slice[j]); */
  ////////////////////////////////////////////////////////////////////////////
  //Get the intervals with respect to slices. Total number of t_slice_min(max) - 2 * n_otu - 3. Excluding tips.
  i = 0;
  For(j, 2 * n_otu - 3)
    {
      t_slice_min[j] = t_slice[i];
      t_slice_max[j] = t_slice[i + 1];
      i++;
    }

  /* For(j, 2 * n_otu - 3) printf("\n. The interval number [%d] min [%f] max[%f] \n", j, t_slice_min[j], t_slice_max[j]); */

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
 
  for(j = 1; j <  2 * n_otu - 3; j++)
    {
      indic[(n_otu - 2) * (2 * n_otu - 3) + j] = 0;
    }

  /* For(i, n_otu - 2) */
  /*   { */
  /*     indic[i * (2 * n_otu - 3)] = 0; */
  /*   } */

  /* For(i, n_otu - 1) */
  /*   { */
  /*     indic[i * (2 * n_otu - 3) + 1] = 0; */
  /*   } */
 
  /* printf("\n"); */
  /* For(i, n_otu - 1) */
  /*   { */
  /*     printf(" ['%d'] ", i + n_otu); */
  /*     For(j, 2 * n_otu - 3) */
  /*       { */
  /*         printf(". '%d' ", indic[i * (2 * n_otu - 3) + j]); */
  /*       } */
  /*     printf("\n"); */
  /*   } */


  ////////////////////////////////////////////////////////////////////////////
  //Running through all of the combinations of slices
  int l, tot_num_comb, *cur_slices, *cur_slices_shr, *cur_slices_cpy, *slices_start_node, shr_num_slices;
  phydbl P, *t_cur_slice_min, *t_cur_slice_max;
  phydbl k_part;
  phydbl num, denom, lmbd;
          
  lmbd = tree -> rates -> birth_rate;

  tot_num_comb = 1;
  P = 0.0; 

 
  t_cur_slice_min    = (phydbl *)mCalloc(n_otu - 1, sizeof(phydbl));
  t_cur_slice_max    = (phydbl *)mCalloc(n_otu - 1, sizeof(phydbl));
  cur_slices         = (int *)mCalloc(n_otu - 1, sizeof(int)); //the vector of the current slices with repetition.
  cur_slices_cpy     = (int *)mCalloc(n_otu - 1, sizeof(int));
  slices_start_node  = (int *)mCalloc(n_otu - 1, sizeof(int)); 
  cur_slices_shr     = (int *)mCalloc(n_otu - 1, sizeof(int)); //the vector of the current slices without repetition.



  ////////////////////////////////////////////////////////////////////////////
  //Get the number of slices that can be applied for each node and the vectors of slice numbers for each node. 
  For(i, n_otu - 1)
    {
      k = 0;
      For(j, 2 * n_otu - 3)
        {
          if(indic[i * (2 * n_otu - 3) + j] == 1) 
            {
              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              g_i_node[i * (2 * n_otu - 3) + j] =  (EXP(lmbd * t_slice_max[j]) - EXP(lmbd * t_slice_min[j])) / 
                (EXP(lmbd * t_prior_max[i + n_otu]) - EXP(lmbd * t_prior_min[i + n_otu])); /* printf("\n. g '%f' ", g_i_node[i * (2 * n_otu - 3) + j]); */
              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              slice_numbers[i * (2 * n_otu - 3) + k] = j; /* printf("\n. Node [%d] slice'%d' ", i + n_otu, slice_numbers[i * (2 * n_otu - 3) + k]); */
              n_slice[i]++; 
              k++;
            }
        }
      /* printf(" Number of slices'%d' \n", n_slice[i]); */
    }
  /* printf("\n"); */
  /* For(i, n_otu - 1) */
  /*   { */
  /*     printf(" ['%d'] ", i + n_otu); */
  /*     For(j, 2 * n_otu - 3) */
  /*       { */
  /*         printf(". '%f' ", g_i_node[i * (2 * n_otu - 3) + j]); */
  /*       } */
  /*     printf("\n"); */
  /*   } */

  /* printf("\n"); */
  /* For(i, n_otu - 1) */
  /*   { */
  /*     printf(" ['%d'] ", i + n_otu); */
  /*     For(j, n_slice[i]) */
  /*       { */
  /*         printf(". '%d' ", slice_numbers[i * (2 * n_otu - 3) + j]); */
  /*       } */
  /*     printf("\n"); */
  /*   } */
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////CALCULATING THE MAXIMUM VALUE OF m_i * g_i FOR ONE COMBINATION OF SLICES:////////////////////////////////////////////////////////////////
  /*   int q; */
  /* For(q, 3){ */
  /* printf("\n"); */
  int r;
  phydbl max, K_total;
  phydbl *t_slice_min_f, *t_slice_max_f;
 
  t_slice_min_f    = (phydbl *)mCalloc(n_otu - 1, sizeof(phydbl));
  t_slice_max_f    = (phydbl *)mCalloc(n_otu - 1, sizeof(phydbl));
  For(i, n_otu - 1)
    {
      r = rand()%(n_slice[i]);
      /* printf("\n. node nmb [%d] r = %d \n", i + n_otu, slice_numbers[i * (2 * n_otu - 3) + r]); */
      t_slice_min_f[i] = t_slice_min[slice_numbers[i * (2 * n_otu - 3) + r]];
      t_slice_max_f[i] = t_slice_max[slice_numbers[i * (2 * n_otu - 3) + r]];
      cur_slices[i] = slice_numbers[i * (2 * n_otu - 3) + r];
      /* max = 0.0; */
      /* /\* printf(" ['%d'] ", i + n_otu);  /\\* printf(" Max value '%f' ", max); *\\/ *\/ */
      /* /\* printf("\n. Slice number '%d' Min '%f' Max '%f' \n", cur_slices[i], t_slice_min_f[i], t_slice_max_f[i]); *\/ */
      /* /\* printf("\n"); *\/ */
      /* For(j, 2 * n_otu - 3) */
      /*   { */
      /*     if(max < g_i_node[i * (2 * n_otu - 3) + j]) */
      /*       { */
      /*         t_slice_min_f[i] = t_slice_min[j]; */
      /*         t_slice_max_f[i] = t_slice_max[j]; */
      /*         /\* printf("\n. [1] Min '%f' Max '%f' \n", t_slice_min_f[i], t_slice_max_f[i]); *\/ */
      /*         max = g_i_node[i * (2 * n_otu - 3) + j]; */
      /*         /\* printf("\n. Max value '%f' \n", max); *\/ */
      /*         /\* printf(". '%f' ", g_i_node[i * (2 * n_otu - 3) + j]); *\/ */
      /*         cur_slices[i] = j; */
      /*       } */
      /*   } */
      /* printf(" Max value '%f' ", max); */
      /* printf(" Slice number '%d' Min '%f' Max '%f' \n", cur_slices[i], t_slice_min_f[i], t_slice_max_f[i]); */
      /* printf("\n"); */
    }
  For(i, n_otu - 1)
    {
      cur_slices_cpy[i]    = cur_slices[i];
      slices_start_node[i] = cur_slices[i];
    }

  For(i, n_otu - 1)
    {
      for(j = i + 1; j < n_otu - 1; j++)
        {
          if(cur_slices[i] == cur_slices[j]) cur_slices[j] = -1;
        }
    }
  //For(i, n_otu - 1) printf(" Slice number'%d' \n", cur_slices[i]);
  
  shr_num_slices = 0;
  For(i, n_otu -1)
    {
      if(cur_slices[i] >= 0)
        {
          cur_slices_shr[shr_num_slices] = cur_slices[i];
          shr_num_slices++;
        }
    }
  /* printf("\n");  */
  /* For(i, shr_num_slices) printf("\n. Slice number'%d' \n", cur_slices_shr[i]);  */
  /* printf("\n");  */
  int result_1;
  int *root_nodes;
  int c = 0;
  /* phydbl n1, n2; */

  root_nodes = (int *)mCalloc(n_otu - 1, sizeof(int));
  
  result_1 = TRUE;
  
  Check_Time_Slices(tree -> n_root, tree -> n_root -> v[1], &result_1, t_slice_min_f, t_slice_max_f, tree);
  Check_Time_Slices(tree -> n_root, tree -> n_root -> v[2], &result_1, t_slice_min_f, t_slice_max_f, tree);
  
  /* printf("\n. [START] Result '%d' \n", result_1); */

  if(result_1) c++; 

  K_total = 0.0;
  
  if(result_1 != TRUE) K_total = 0.0;
  else
    {
      int n_1, n_2;
      /* int *root_nodes; */
      int num_elem;
      
      num_elem = 0;
      
      /* root_nodes = (int *)mCalloc(n_otu - 1, sizeof(int)); */
      
      For(i, shr_num_slices)
        {
          Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[1], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_slice_min_f, t_slice_max_f, tree);
          Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[2], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_slice_min_f, t_slice_max_f, tree);   /* printf(" Max value '%f' ", max); */
      /* printf("\n. Slice number '%d' Min '%f' Max '%f' \n", cur_slices[i], t_slice_min_f[i], t_slice_max_f[i]); */
      /* printf("\n"); */
        }
      /* for(i = 0; i < n_otu - 2; i++) printf("\n. [START] Node [%d] Slice_min [%f] Slice_max [%f] \n", i + n_otu, t_slice_min_f[i], t_slice_max_f[i]); */
      For(j, num_elem)
        {
          n_1 = 0;
          n_2 = 0;
          /* n1 = 0.0; */
          /* n2 = 0.0; */
          Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[1], &n_1, t_slice_min_f, t_slice_max_f, tree);
          Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[2], &n_2, t_slice_min_f, t_slice_max_f, tree);
          /* printf("\n. n_1 [%d] n_2 [%d]\n", n_1, n_2); */
          /* printf("\n. n_1+n_2 [%f]  n_1+n_2+1 [%f]  n_1 [%f]  n_2 [%f]\n", LOG(Factorial(n_1 + n_2)), LOG(Factorial(n_1 + n_2 + 1)), LOG(Factorial(n_1)), LOG(Factorial(n_2))); */
          /* K_total = K_total + LOG(Factorial(n_1 + n_2)) - LOG(n_1 + n_2 + 1) - LOG(Factorial(n_1 + n_2)) - LOG(Factorial(n_1)) - LOG(Factorial(n_2)); */
          /* K_total = K_total * (phydbl)Factorial(n_1 + n_2) / ((phydbl)Factorial(n_1 + n_2 + 1) * (phydbl)Factorial(n_1) * Factorial(n_2)); */
          /* K_total = K_total + LOG(1) - LOG(n_1 + n_2 + 1) - LOG(Factorial(n_1)) - LOG(Factorial(n_2)); */
          K_total = K_total + LOG(1) - LOG(n_1 + n_2 + 1) - LnGamma(n_1 + 1) - LnGamma(n_2 + 1);
          /* for(i = 1; i < n_1 +1; i++) n1 = n1 + LOG(i); */
          /* for(i = 1; i < n_2 +1; i++) n2 = n2 + LOG(i); */
          /* K_total = K_total + LOG(1) - LOG(n_1 + n_2 + 1) - n1 - n2; */
          /* printf("\n. K_total [%f] \n", K_total); */
          /* printf("\n. [START] LOG(m_i) [%f] \n", LOG(1) - LOG(n_1 + n_2 + 1) - LnGamma(n_1 + 1) - LnGamma(n_2 + 1)); */
        }
      /* printf("\n. [START] LOG(m_i) [%f] \n", K_total); */
      /* printf("\n. K_total [%f] \n", K_total); */
 
 
      /* num = 1; */
      /* denom = 1; */
      num = 0.0;
      denom = 0.0;
      /* lmbd = 4.0; */

 

      for(i = n_otu; i < 2 * n_otu - 2; i++) if(Are_Equal(t_prior_min[tree -> n_root -> num], t_prior_min[i], 1.E-10)) t_prior_min[i] = tree -> rates -> nd_t[tree -> n_root -> num];
      For(j, n_otu - 2)
        {
          num = num + LOG(EXP(lmbd * t_slice_max_f[j]) - EXP(lmbd * t_slice_min_f[j]));
          /* printf("\n. num [%f] \n", num); */
        }
       for(j = n_otu; j < 2 * n_otu - 2; j++)
        {
          denom = denom + LOG(EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j]));
          /* printf("\n. denom [%f] \n", denom); */
        }

       /* printf("\n. [START] LOG(g_i) [%f] \n", num - denom); */
       /* K_total = (K_total * num) / denom; */
       K_total = EXP(K_total + num - denom);
       /* printf("\n. [START] sum(m_i * g_i) = m_i * g_i [%f] \n", K_total); */
       if(isinf(K_total) || isnan(K_total)) 
         {
           printf("\n. [1] LMBD %f \n", lmbd);
           For(j, n_otu - 2) printf("\n. [1] EXPdif %f LOG(ESPdif) %f t_slice_min %f t_slice_max %f \n", EXP(lmbd * t_slice_max_f[j]) - EXP(lmbd * t_slice_min_f[j]), LOG(EXP(lmbd * t_slice_max_f[j]) - EXP(lmbd * t_slice_min_f[j])), t_slice_min_f[j], t_slice_max_f[j]);
           for(j = n_otu; j < 2 * n_otu - 2; j++)  printf("\n. [2] EXPdif %f LOG(ESPdif) %f t_prior_min %f t_prior_max %f \n", EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j]), LOG(EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j])), t_prior_min[j], t_prior_max[j]);
           PhyML_Printf("\n. K_total=%f \n", K_total);
           PhyML_Printf("\n. n_1=%d n_2=%d \n", n_1, n_2);
           PhyML_Printf("\n. num=%f denom=%f \n", num, denom);
           PhyML_Printf("\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
           Warn_And_Exit("\n");
         }
       /* printf("\n. K of the tree for starting combination of slices [%f] \n", K_total); */
    }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /* printf("\n. ____________________________________________________________________________________________ \n"); */
  phydbl K_part, K_max;
  int m, n, count; /* count = 0; */

  K_max = K_total;
  /* printf("\n. K_total [%f] \n", K_total); */
  do
    {
      count = 0;
      For(m, n_otu - 1)
        {
          /* printf("\n. NODE NUMBER ['%d'] \n", m + n_otu); */
          For(n, 2 * n_otu - 3)
            {
              /* printf("\n"); */
              /* K_part = 1; */
              K_part = 0.0;
              if(g_i_node[m * (2 * n_otu - 3) + n] > 0)
                {
                  t_slice_min_f[m] = t_slice_min[n];
                  t_slice_max_f[m] = t_slice_max[n];
                  cur_slices_cpy[m] = n;
                  /* For(i, n_otu - 1) printf(" [%d] ", cur_slices_cpy[i]); */
                  /* printf("\n"); */
                  For(i, n_otu - 1) if(!Are_Equal(cur_slices_cpy[i], slices_start_node[i], 1.E-10))
                    {
                      /* For(i, n_otu - 1) printf(" [%d] ", cur_slices_cpy[i]); */
                      /* printf("\n"); */
                      /* For(i, n_otu - 2) printf("\n. [2] Node number [%d] Min '%f' Max '%f' \n", i + n_otu, t_slice_min_f[i], t_slice_max_f[i]); */
                      For(i, n_otu - 1) cur_slices[i] = cur_slices_cpy[i];
                      /* For(i, n_otu) printf("\n. [1] Slice numbers [%d] \n", cur_slices[i]); */
                      For(i, n_otu - 1)
                        {
                          for(j = i + 1; j < n_otu - 1; j++)
                            {
                              if(cur_slices[i] == cur_slices[j]) cur_slices[j] = -1;
                            }
                        }
                      
                      shr_num_slices = 0;
                      For(i, n_otu -1)
                        {
                          if(cur_slices[i] >= 0)
                            {
                              cur_slices_shr[shr_num_slices] = cur_slices[i];
                              shr_num_slices++;
                            }
                        }
                      /* For(i, shr_num_slices) printf("\n. [2] Slice numbers [%d] \n", cur_slices_shr[i]); */
                      
                      result_1 = TRUE;
                      
                      Check_Time_Slices(tree -> n_root, tree -> n_root -> v[1], &result_1, t_slice_min_f, t_slice_max_f, tree);
                      Check_Time_Slices(tree -> n_root, tree -> n_root -> v[2], &result_1, t_slice_min_f, t_slice_max_f, tree);
                      
                      /* printf("\n. [CONT] Result '%d' \n", result_1); */

                      if(result_1) c++; 
                      
                      if(result_1 != TRUE) K_part = 0.0;
                      else
                        {
                          int n_1, n_2;
                          /* int *root_nodes; */
                          int num_elem;
                          
                          num_elem = 0;
                          
                          /* root_nodes = (int *)mCalloc(n_otu - 1, sizeof(int)); */
                          
                          For(i, shr_num_slices)
                            {
                              Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[1], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_slice_min_f, t_slice_max_f, tree);
                              Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[2], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_slice_min_f, t_slice_max_f, tree);
                            }
                          /* for(i = 0; i < n_otu - 2; i++) printf("\n. [CONT] Node [%d] Slice_min [%f] Slice_max [%f] \n", i + n_otu, t_slice_min_f[i], t_slice_max_f[i]); */
                          For(j, num_elem)
                            {
                              n_1 = 0;
                              n_2 = 0;
                              /* n1 = 0.0; */
                              /* n2 = 0.0; */
                              Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[1], &n_1, t_slice_min_f, t_slice_max_f, tree);
                              Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[2], &n_2, t_slice_min_f, t_slice_max_f, tree);
                              /* printf("\n. n_1 [%d] n_2 [%d]\n", n_1, n_2); */
                              /* K_part = K_part * Factorial(n_1 + n_2) / ((phydbl)Factorial(n_1 + n_2 + 1) * Factorial(n_1) * Factorial(n_2)); */
                              /* K_part = K_part + LOG(Factorial(n_1 + n_2)) - LOG(n_1 + n_2 + 1) - LOG(Factorial(n_1 + n_2)) - LOG(Factorial(n_1)) - LOG(Factorial(n_2)); */
                              /* K_part = K_part + LOG(1) - LOG(n_1 + n_2 + 1) - LOG(Factorial(n_1)) - LOG(Factorial(n_2)); */
                              K_part = K_part + LOG(1) - LOG(n_1 + n_2 + 1) - LnGamma(n_1 + 1) - LnGamma(n_2 + 1);
                              /* for(i = 1; i < n_1 +1; i++) n1 = n1 + LOG(i); */
                              /* for(i = 1; i < n_2 +1; i++) n2 = n2 + LOG(i); */
                              /* K_total = K_total + LOG(1) - LOG(n_1 + n_2 + 1) - n1 - n2; */
                              /* K_total = K_total + LOG(1) - LOG(n_1 + n_2 + 1) - LnGamma(n_1 + 1) - LnGamma(n_2 + 1); */
                              /* printf("\n. [CONT] LOG(m_i) [%f] \n", LOG(1) - LOG(n_1 + n_2 + 1) - LnGamma(n_1 + 1) - LnGamma(n_2 + 1)); */
                            }
                          /* printf("\n. [CONT] LOG(m_i) [%f] \n", K_part); */
                          /* printf("\n. K_part [%f] \n", K_part); */
                          
                          num = 0.0;
                          denom = 0.0;
                          /* num = 1; */
                          /* denom = 1; */
                          /* lmbd = 4.0; */

                        
                           
                          for(i = n_otu; i < 2 * n_otu - 2; i++) if(Are_Equal(t_prior_min[tree -> n_root -> num], t_prior_min[i], 1.E-10)) t_prior_min[i] = tree -> rates -> nd_t[tree -> n_root -> num];
                          For(j, n_otu - 2)
                            {
                              num = num + LOG(EXP(lmbd * t_slice_max_f[j]) - EXP(lmbd * t_slice_min_f[j]));
                            }
                          for(j = n_otu; j < 2 * n_otu - 2; j++)
                            {
                              denom = denom + LOG(EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j]));
                            }

                          /* printf("\n. [CONT] LOG(g_i) [%f] \n", num - denom); */
                          /* K_part = (K_part * num) / denom; */                   

                          K_part = EXP(K_part + num - denom);

                          /* printf("\n. [START] m_i * g_i [%f] \n", K_part); */
                          if(isinf(K_part) || isnan(K_part)) 
                            {
                              printf("\n. [2] LMBD %f \n", lmbd);
                              For(j, n_otu - 2) printf("\n. [1] EXPdif %f LOG(ESPdif) %f t_slice_min %f t_slice_max %f \n", EXP(lmbd * t_slice_max_f[j]) - EXP(lmbd * t_slice_min_f[j]), LOG(EXP(lmbd * t_slice_max_f[j]) - EXP(lmbd * t_slice_min_f[j])), t_slice_min_f[j], t_slice_max_f[j]);
                              for(j = n_otu; j < 2 * n_otu - 2; j++)  printf("\n. [2] EXPdif %f LOG(ESPdif) %f t_prior_min %f t_prior_max %f \n", EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j]), LOG(EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j])), t_prior_min[j], t_prior_max[j]);
                              PhyML_Printf("\n. K_part=%f \n", K_part);
                              PhyML_Printf("\n. n_1=%d n_2=%d \n", n_1, n_2);
                              PhyML_Printf("\n. num=%f denom=%f \n", num, denom);
                              PhyML_Printf("\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
                              Warn_And_Exit("\n");
                            }
                          K_total = K_total + K_part;
                          /* printf("\n. K_max [%f] K_part [%f] \n", K_max, K_part); */
                          /* printf("\n. [CONT] sum(m_i * g_i) [%f] \n", K_total); */
                          if(K_max < K_part)
                            {
                              K_max = K_part;
                              /* For(i, n_otu - 1) slices_start_node[i] = cur_slices_cpy[i]; */
                              count++;
                            }
                          /* printf("\n. K_part [%f] \n", K_part); */
                          /* printf("\n. K_total [%f] \n", K_total); */
                          /* printf("\n. [1] Normolizing constant [%f] \n", 1 / K_total); */
                          /* count++; */
                          /* printf("\n"); */
                          /* printf("\n. Count [%d] \n", count); */
                          /* printf("\n"); */
                        }
                    }
                }
            }
          /* For(i, n_otu - 1) slices_start_node[i] = cur_slices_cpy[i]; */
          /* printf("\n. [2] Normolizing constant [%f] \n", 1 / (K_total)); */
          /* printf("\n. ____________________________________________________________________________________________ \n");       */
          /* Exit("\n"); */
          For(i, n_otu - 1) slices_start_node[i] = cur_slices_cpy[i];
        }
      if(Are_Equal(count, 0, 1.E-10)) break;
    }
  while(1);
  
  /* printf("\n. [%d] Normolizing constant [%f] \n", q+1, 1 / (K_total)); */
  /* printf("\n. [APPROX] Total number of TRUE comb of slices used for approximation [%d] \n", c); */
  /* printf("\n. [APPROX] Approximated constant [%f] \n", 1 / (K_total)); */
  /* Exit("\n"); */
  /* printf("\n. ____________________________________________________________________________________________ \n"); */
  /* } */
  /* n1 = 6; */
  /* printf("\n. LogOfFact %f \n", LOG(Factorial(n1))); */
  /* printf("\n. LnGamma %f \n", LnGamma(n1+1)); */
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /* For(i, n_otu - 1) */
  /*   { */
  /*     /\* printf("\n. [2] Number of slices'%d' \n", n_slice[i]); *\/ */
  /*     tot_num_comb = tot_num_comb * n_slice[i]; */
  /*     /\* printf("\n. [2] Total number of combinations of slices [%d] \n", tot_num_comb); *\/ */
  /*   } */

  /* /\* printf("\n. Total number of combinations of slices [%d] \n", tot_num_comb); *\/ */
  /* int p = 0; */
  /* For(k, tot_num_comb) */
  /*   { */
  /*     shr_num_slices = 0; */
  /*     /\* printf("\n"); *\/ */
  /*     For(i, n_otu - 1) //node number i + n_otu */
  /*       { */
  /*         /\* printf(" [%d] ", i + n_otu); *\/ */
  /*         l = (k % Number_Of_Comb_Slices(i, n_otu - 1, n_slice)) / Number_Of_Comb_Slices(i+1, n_otu - 1, n_slice); //printf(" Slice number'%d' ", slice_numbers[i * (2 * n_otu - 3) + l]); */
  /*         t_cur_slice_min[i] = t_slice_min[slice_numbers[i * (2 * n_otu - 3) + l]]; /\* printf(" '%f' ", t_cur_slice_min[i]); *\/ */
  /*         t_cur_slice_max[i] = t_slice_max[slice_numbers[i * (2 * n_otu - 3) + l]]; /\* printf(" '%f' ", t_cur_slice_max[i]); *\/ */
  /*         cur_slices[i] = slice_numbers[i * (2 * n_otu - 3) + l]; */
  /*         /\* printf("\n"); *\/ */
  /*       } */
  /*     //printf("\n"); */
  /*     //For(i, n_otu - 1) printf(" Slice number'%d' ", cur_slices[i]); */
  /*     //printf("\n"); */

  /*     /////////////////////////////////////////////////////////////////////////// */
  /*     //Taking away duplicated slices */
  /*     For(i, n_otu - 1) */
  /*       { */
  /*         for(j = i + 1; j < n_otu - 1; j++) */
  /*           { */
  /*             if(cur_slices[i] == cur_slices[j]) cur_slices[j] = -1; */
  /*           } */
  /*       } */
  /*     //For(i, n_otu - 1) printf(" Slice number'%d' \n", cur_slices[i]); */

  /*     /////////////////////////////////////////////////////////////////////////// */
  /*     //Getting a vector of all of the slices without duplicates. */
  /*     For(i, n_otu -1) */
  /*       { */
  /*         if(cur_slices[i] >= 0) */
  /*           { */
  /*             cur_slices_shr[shr_num_slices] = cur_slices[i]; */
  /*             shr_num_slices++; */
  /*           } */
  /*       } */
  /*     //printf("\n"); */
  /*     //For(i, shr_num_slices) printf("\n. Slice number'%d' \n", cur_slices_shr[i]); */
  /*     //printf("\n"); */

  /*     //////////////////////////////////////////////////////////////////////////// */
  /*     //Check for the time slices to be set properly */
  /*     int result; */

  /*     result = TRUE; */

  /*     Check_Time_Slices(tree -> n_root, tree -> n_root -> v[1], &result, t_cur_slice_min, t_cur_slice_max, tree); */
  /*     Check_Time_Slices(tree -> n_root, tree -> n_root -> v[2], &result, t_cur_slice_min, t_cur_slice_max, tree); */
  /*     //printf("\n. '%d' \n", result); */

  /*     if(result) p++; */

  /*     //For(i, n_otu - 1) printf("\n. Node [%d] min [%f] max [%f] \n", i + n_otu, t_cur_slice_min[i], t_cur_slice_max[i]); */

  /*     //////////////////////////////////////////////////////////////////////////// */
  /*     //Calculating k_part */

  /*     k_part = 1.0; */

  /*     if(result != TRUE) k_part = 0.0; */
  /*     else */
  /*       { */
  /*         int n_1, n_2; */

  /*         //////////////////////////////////////////////////////////////////////////// */
  /*         //Getting the root node in a given slice */
  /*         int *root_nodes; */
  /*         int num_elem; */
           
  /*         num_elem = 0; */
           
  /*         root_nodes = (int *)mCalloc(n_otu - 1, sizeof(int)); */
          
  /*         //printf("\n. Number of slices shrinked [%d] \n", shr_num_slices); */
  /*         /\* for(i = 0; i < n_otu - 2; i++) printf("\n. [EXACT] Node [%d] Slice_min [%f] Slice_max [%f] \n", i + n_otu, t_cur_slice_min[i], t_cur_slice_max[i]);     *\/    */
  /*         For(i, shr_num_slices) */
  /*           { */
  /*             //printf("\n. The number of the shrinked interval [%d] min [%f] max [%f] \n", cur_slices_shr[i], t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]]); */
  /*             Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[1], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_cur_slice_min, t_cur_slice_max, tree); */
  /*             Search_Root_Node_In_Slice(tree -> n_root, tree -> n_root -> v[2], root_nodes, &num_elem, t_slice_min[cur_slices_shr[i]], t_slice_max[cur_slices_shr[i]], t_cur_slice_min, t_cur_slice_max, tree); */
  /*           } */
  /*         //printf("\n. Number of elements in a vector of the nodes [%d] \n", num_elem); */
  /*         //For(j, num_elem) printf("\n. Root node number [%d] \n", root_nodes[j] + n_otu); */
  /*         For(j, num_elem) */
  /*           { */
  /*             //printf("\n. Root node number [%d] \n", tree -> a_nodes[root_nodes[j] + n_otu] -> num); */
              
  /*             n_1 = 0; */
  /*             n_2 = 0; */
  /*             Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[1], &n_1, t_cur_slice_min, t_cur_slice_max, tree); */
  /*             Number_Of_Nodes_In_Slice(tree -> a_nodes[root_nodes[j] + n_otu], tree -> a_nodes[root_nodes[j] + n_otu] -> v[2], &n_2, t_cur_slice_min, t_cur_slice_max, tree); */
  /*             //printf("\n. n_1 [%d] n_2 [%d]\n", n_1, n_2); */
  /*             k_part = k_part * Factorial(n_1 + n_2) / ((phydbl)Factorial(n_1 + n_2 + 1) * Factorial(n_1) * Factorial(n_2)); */
  /*           } */

  /*         /\* printf("\n. [EXACT] m_i [%f] \n", k_part);    *\/ */
  /*         /\* printf("\n. k_part [%f] \n", k_part);      *\/ */
  /*         //////////////////////////////////////////////////////////////////////////// */
  /*         //Calculating PRODUCT over all of the time slices k_part * (exp(-lmbd*l) - exp(-lmbd*u))/(exp(-lmbd*l) - exp(-lmbd*u)) */
         
  /*         /\* lmbd = tree -> rates -> birth_rate; *\/ */
  /*         num = 1; */
  /*         denom = 1; */
  /*         //lmbd = 4.0; */
  /*         for(i = n_otu; i < 2 * n_otu - 2; i++) if(Are_Equal(t_prior_min[tree -> n_root -> num], t_prior_min[i], 1.E-10)) t_prior_min[i] = tree -> rates -> nd_t[tree -> n_root -> num]; */
  /*         For(j, n_otu - 2) */
  /*           { */
  /*             num = num * (EXP(lmbd * t_cur_slice_max[j]) - EXP(lmbd * t_cur_slice_min[j])); */
  /*           } */
  /*         for(j = n_otu; j < 2 * n_otu - 2; j++) */
  /*           { */
  /*             denom = denom * (EXP(lmbd * t_prior_max[j]) - EXP(lmbd * t_prior_min[j])); */
  /*           } */
  /*         /\* printf("\n. [EXACT] g_i [%f] \n", num / denom);    *\/ */
  /*         k_part = (k_part * num) / denom; */
  /*         /\* printf("\n. [EXACT] m_i * g_i [%f] \n", k_part);  *\/ */
  /*         /\* printf("\n. [EXACT] m_i * g_i [%f] comb numb [%d] \n", k_part, k);   *\/ */
  /*         /\* if(k_part > 0.0001) for(i = 0; i < n_otu - 2; i++) printf("\n. [EXACT] Node [%d] Slice_min [%f] Slice_max [%f] \n", i + n_otu, t_cur_slice_min[i], t_cur_slice_max[i]);      *\/ */
  /*        } */
  /*     P = P + k_part; */
  /*     /\* printf("\n. [EXACT] sum(m_i * g_i) [%f] \n", P);    *\/ */
  /*   } */
  /* //printf("\n. [P] of the tree for one combination of slices [%f] \n", P); */
  /* K = 1 / P; */
  /* printf("\n. [EXACT] Total number of TRUE comb of slices [%d] \n", p); */
  /* printf("\n. [EXACT] Exact constant [%f] \n", K); */
  /* /\* printf("\n. LOG(g_i) [%f] \n",  LOG_g_i(lmbd, -10.0, -20.0, 0.0, -30.0)); *\/ */
  /* /\* int a, b, d; *\/ */
  /* /\* a = 1; *\/ */
  /* /\* b = 3; *\/ */
  /* /\* d = 1; *\/ */

  /* /\* printf("\n. [%d] \n",  (a % b) / d); *\/ */
  /* /\* printf("\n. LOG(g_i) [%f] \n",  LOG_g_i(lmbd, t_slice_max, t_slice_min, t_prior_max, t_prior_min)); *\/ */
  /* /\* Exit("\n"); *\/ */
  /* printf("\n. ____________________________________________________________________________________________ \n"); */
 
  ////////////////////////////////////////////////////////////////////////////
  free(t_cur_slice_min);
  free(t_cur_slice_max);
  free(cur_slices);
  free(cur_slices_cpy);
  free(slices_start_node);
  free(cur_slices_shr); 
  free(t_slice);        
  free(t_slice_min);    
  free(t_slice_max); 
  free(t_slice_min_f);    
  free(t_slice_max_f);     
  free(indic);          
  free(slice_numbers); 
  free(root_nodes);   
  free(n_slice);
  free(g_i_node);    
  //free(tree -> K);
 
  /* return(K); */
  return(1 / K_total);
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
          if(calib -> all_applies_to[k] -> num)
            {
              t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower);
              t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper);
              t_has_prior[calib -> all_applies_to[k] -> num] = YES;
            }
          /* else */
          /*   { */
          /*     if(calib -> next) calib = calib -> next; */
          /*     else break; */
          /*   } */
          if(calib -> next) calib = calib -> next;
          else break;
        }
      while(calib);
      int result;
      TIMES_Set_All_Node_Priors_S(&result, tree);
      //for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) printf("\n. [1] Node [%d] min [%f] max[%f]\n", j, tree -> rates -> t_prior_min[j], tree -> rates -> t_prior_max[j]);
      //tree -> rates -> birth_rate = 4.0;
      K[i] = Slicing_Calibrations(tree);     
      /* PhyML_Printf("\n. Number [%d] normolizing constant [%f] \n", i+1, K[i]); */
      while(calib -> prev) calib = calib -> prev;
    }
  return(K);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Sets a vector of the partial probabilities for each combination of calibrations
void TIMES_Calib_Partial_Proba(t_tree *tree)
{

  phydbl *times_partial_proba, proba, *t_prior_min, *t_prior_max; 
  int i, j, k, tot_num_comb;
  t_cal *calib;
  short int *t_has_prior;

  proba = 0.0;
 
  times_partial_proba = tree -> rates -> times_partial_proba; 
  calib = tree -> rates -> calib;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;

  tot_num_comb = Number_Of_Comb(calib);

  For(i, tot_num_comb)
    { 
      /* printf("\n\n"); */
      times_partial_proba[i] = 1.0;
      for(j = tree -> n_otu; j < 2 * tree -> n_otu - 1; j++) 
        {
          t_prior_min[j] = -BIG;
          t_prior_max[j] = BIG;
          t_has_prior[j] = NO; 
        }
      do
        {
          k = (i % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next);
          if(calib -> all_applies_to[k] -> num)
            {
              /* printf("\n. Node number [%d] \n", calib -> all_applies_to[k] -> num); */
              t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower); /* printf("\n. [1] Min [%f] \n",  t_prior_min[calib -> all_applies_to[k] -> num]); */
              t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper); /* printf("\n. [1] Max [%f] \n",  t_prior_max[calib -> all_applies_to[k] -> num]); */
              t_has_prior[calib -> all_applies_to[k] -> num] = YES;
              proba = calib -> proba[calib -> all_applies_to[k] -> num]; 
              times_partial_proba[i] *= proba;
              /* printf("\n. [1] Proba [%f] \n", proba); */
              /* if((t_prior_min[calib -> all_applies_to[k] -> num] > t_prior_max[calib -> all_applies_to[k] -> num])) times_partial_proba[i] = 0.0;  */
              /* else times_partial_proba[i] *= calib -> proba[calib -> all_applies_to[k] -> num];  */
            }
          else
            {
              /* printf("\n. k [%d] \n", k);    */         
              proba = calib -> proba[2 * tree -> n_otu - 1]; 
              times_partial_proba[i] *= proba;
              /* printf("\n. [2] Proba [%f] \n", proba); */
              /* if(calib -> next) calib = calib -> next; */
              /* else break; */
            }
          
          if(calib -> next) calib = calib -> next;
          else break;
        }
      while(calib);

      int result;

      result = TRUE;

      /* printf("\n. [3] Partial Proba [%f] \n", times_partial_proba[i]); */

      TIMES_Set_All_Node_Priors_S(&result, tree);

      if(result != TRUE) times_partial_proba[i] = 0; /* printf("\n. [4] Partial Proba [%f] \n", times_partial_proba[i]); */
      /* times_partial_proba[i] = 1.0; */
      /* do */
      /*   { */
      /*     k = (i % Number_Of_Comb(calib)) / Number_Of_Comb(calib -> next); */
      /*     if(calib -> all_applies_to[k] -> num) proba = calib -> proba[calib -> all_applies_to[k] -> num];  */
      /*     else proba = calib -> proba[2 * tree -> n_otu - 1]; */
      /*     times_partial_proba[i] *= proba; */
      /*     if(calib -> next) calib = calib -> next; */
      /*     else break; */
      /*   } */
      /* while(calib); */
      while(calib -> prev) calib = calib -> prev;
      /* printf("\n\n"); */
    }

  phydbl sum_proba;
  sum_proba = 0.0;
  /* For(i, tot_num_comb) printf("\n. [1] Partial Proba [%f] \n", times_partial_proba[i]);  */
  For(i, tot_num_comb) sum_proba += times_partial_proba[i];
  if(!Are_Equal(sum_proba, 1.0, 1.E-10)) 
    {
      For(i, tot_num_comb) times_partial_proba[i] = times_partial_proba[i] / sum_proba;
    } 
  /* For(i, tot_num_comb) printf("\n. [2] Partial Proba [%f] \n", times_partial_proba[i]);  */
  /* Exit("\n"); */
}



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
          return;
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
//Function calculates the TOTAL number of calibration combinations, 
//given the number of nodes to which each calibartion applies to.
int Number_Of_Calib(t_cal *calib)
{

  int num_calib;


  num_calib = 0;
  do
    {
      num_calib++;
      if(calib -> next) calib = calib -> next;
      else break;
    }
  while(calib);
  return(num_calib);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//Function sets current calibartion in the following way:
//Suppose we have a vector of calibrations C=(C1, C2, C3), each calibration  
//applies to a set of nodes. we can reach each node number through the indeces (corresponds 
//to the number the information was read). C1={0,1,2}, C2={0,1}, C3={0};
//The total number of combinations is 3*2*1=6. The first combination with row number 0 
//will be {0,0,0}, the second row will be {0,1,0} and so on. Calling the node numbers with 
//the above indeces will return current calibration. Also sets the vector of the probabilities
//for current calibration combination.   
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
      if(calib -> all_applies_to[k] -> num) 
        {
          t_prior_min[calib -> all_applies_to[k] -> num] = MAX(t_prior_min[calib -> all_applies_to[k] -> num], calib -> lower);
          t_prior_max[calib -> all_applies_to[k] -> num] = MIN(t_prior_max[calib -> all_applies_to[k] -> num], calib -> upper);     
          t_has_prior[calib -> all_applies_to[k] -> num] = YES; 
          curr_nd_for_cal[i] = calib -> all_applies_to[k] -> num;
          i++;
        }
      else 
        {
          if(calib->next) calib = calib->next;
          else break;    
        }     
      if(calib->next) calib = calib->next;
      else break;    
    }
  while(calib);
  //while(calib -> prev) calib = calib -> prev;
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
  phydbl *t_prior_min, *t_prior_max; //*times_partial_proba;
  short int *t_has_prior;
  t_cal *cal;

  tot_num_cal = tree -> rates -> tot_num_cal;
  t_prior_min = tree -> rates -> t_prior_min;
  t_prior_max = tree -> rates -> t_prior_max;
  t_has_prior = tree -> rates -> t_has_prior;
  //times_partial_proba = tree -> rates -> times_partial_proba;
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
      t_prior_min[j] = -BIG;
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


  t_low = MAX(t_prior_min[d -> num], nd_t[a -> num]);
  t_up  = t_prior_max[d -> num];  

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
          /* t_prior_min[d -> num] = t_low;  */
          /* t_prior_max[d -> num] = t_up;   */
        }    
      
      For(i,3) 
        if((d -> v[i] != d -> anc) && (d -> b[i] != tree -> e_root)) 
          Update_Times_Down_Tree(d, d -> v[i], L_Hastings_ratio, tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

xml_node *XML_Search_Node_Attribute_Value_Clade(char *attr_name, char *value, int skip, xml_node *node)
{
  xml_node *match;

  //printf("\n. Node name [%s] Attr name [%s] Attr value [%s] \n", node -> name, attr_name, value);
  if(!node)
    {
      PhyML_Printf("\n== node: %p attr: %p",node,node?node->attr:NULL);
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");         
    }
  
  match = NULL;
  if(skip == NO && node -> attr)
    {
      xml_attr *attr;
      attr = node -> attr;
      do
        {
          if(!strcmp(attr -> name, attr_name) && !strcmp(attr -> value, value))
            {
              match = node;
              break;
            }
          attr = attr->next;
          if(!attr) break;
        }
      while(1);
    }
  if(match) return(match);
  if(node -> next && !match)
    {
      match = XML_Search_Node_Attribute_Value_Clade(attr_name, value, NO, node -> next);
      return match;
    }
  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
char **XML_Reading_Clade(xml_node *n_clade, t_tree *tree)
{
  int i, /* j, */ n_otu;
  char **clade;

  i = 0;
  n_otu = tree -> n_otu;

  clade = (char **)mCalloc(n_otu, sizeof(char *));
  /* For(j, n_otu) */
  /*   { */
  /*     clade[j] = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
  /*   } */

  if(n_clade)
    {
      do
        {
          clade[i] =  n_clade -> attr -> value; 
          i++;
          if(n_clade -> next) n_clade = n_clade -> next;
          else break;
        }
      while(n_clade);
    }
  else
    {
      PhyML_Printf("== Clade is empty. \n");
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return(clade);                          
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int XML_Number_Of_Taxa_In_Clade(xml_node *n_clade)
{
  int clade_size = 0;
  if(n_clade)
    {
      do
        {
          clade_size++; 
          if(n_clade -> next) n_clade = n_clade -> next;
          else break;
        }
      while(n_clade);
    }
  else
    {
      PhyML_Printf("==Clade is empty. \n");
      PhyML_Printf("\n== Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }
  return(clade_size);
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void TIMES_Set_All_Node_Priors_S(int *result, t_tree *tree)
{
  int i;
  phydbl min_prior;


  /* Set all t_prior_max values */
  TIMES_Set_All_Node_Priors_Bottom_Up_S(tree->n_root,tree->n_root->v[2], result, tree);
  TIMES_Set_All_Node_Priors_Bottom_Up_S(tree->n_root,tree->n_root->v[1], result, tree);

  tree->rates->t_prior_max[tree->n_root->num] =
    MIN(tree->rates->t_prior_max[tree->n_root->num],
        MIN(tree->rates->t_prior_max[tree->n_root->v[2]->num],
            tree->rates->t_prior_max[tree->n_root->v[1]->num]));


  /* Set all t_prior_min values */
  if(!tree->rates->t_has_prior[tree->n_root->num])
    {
      min_prior = 1.E+10;
      For(i,2*tree->n_otu-2)
        {
          if(tree->rates->t_has_prior[i])
            {
              if(tree->rates->t_prior_min[i] < min_prior)
        	min_prior = tree->rates->t_prior_min[i];
            }
        }
      tree->rates->t_prior_min[tree->n_root->num] = 2.0 * min_prior;
      /* tree->rates->t_prior_min[tree->n_root->num] = 10. * min_prior; */
    }
  
  if(tree->rates->t_prior_min[tree->n_root->num] > 0.0)
    {
      /* PhyML_Printf("\n== Failed to set the lower bound for the root node."); */
      /* PhyML_Printf("\n== Make sure at least one of the calibration interval"); */
      /* PhyML_Printf("\n== provides a lower bound."); */
      /* Exit("\n"); */
      *result = FALSE;
      /* return; */
    }


  TIMES_Set_All_Node_Priors_Top_Down_S(tree->n_root,tree->n_root->v[2], result, tree);
  TIMES_Set_All_Node_Priors_Top_Down_S(tree->n_root,tree->n_root->v[1], result, tree);

  /* Get_Node_Ranks(tree); */
  /* TIMES_Set_Floor(tree); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Set_All_Node_Priors_Bottom_Up_S(t_node *a, t_node *d, int *result, t_tree *tree)
{
  int i;
  phydbl t_sup;

  if(d->tax) return;
  else 
    {
      t_node *v1, *v2; /* the two sons of d */

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      TIMES_Set_All_Node_Priors_Bottom_Up_S(d,d->v[i], result, tree);	      
	    }
	}
      
      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
      
      if(tree->rates->t_has_prior[d->num] == YES)
	{
	  t_sup = MIN(tree->rates->t_prior_max[d->num],
		      MIN(tree->rates->t_prior_max[v1->num],
			  tree->rates->t_prior_max[v2->num]));

	  tree->rates->t_prior_max[d->num] = t_sup;

	  /* if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num] || tree->rates->nd_t[d->num] < tree->rates->t_prior_min[d->num] */
          /*    || tree->rates->nd_t[d->num] > tree->rates->t_prior_max[d->num]) */
	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      /* PhyML_Printf("\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]); */
	      /* PhyML_Printf("\n. Inconsistency in the prior settings detected at t_node %d",d->num); */
	      /* PhyML_Printf("\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__); */
	      /* Warn_And_Exit("\n"); */
              *result = FALSE;
              /* return; */
	    }
	}
      else
	{
	  tree->rates->t_prior_max[d->num] = 
	    MIN(tree->rates->t_prior_max[v1->num],
		tree->rates->t_prior_max[v2->num]);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TIMES_Set_All_Node_Priors_Top_Down_S(t_node *a, t_node *d, int *result, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      
      if(tree->rates->t_has_prior[d->num] == YES)
	{
	  tree->rates->t_prior_min[d->num] = MAX(tree->rates->t_prior_min[d->num],tree->rates->t_prior_min[a->num]);
	  
	  /* if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num]  || tree->rates->nd_t[d->num] < tree->rates->t_prior_min[d->num] */
          /*    || tree->rates->nd_t[d->num] > tree->rates->t_prior_max[d->num]) */
	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      /* PhyML_Printf("\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]); */
	      /* PhyML_Printf("\n. Inconsistency in the prior settings detected at t_node %d",d->num); */
	      /* PhyML_Printf("\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__); */
	      /* Warn_And_Exit("\n"); */
              *result = FALSE;
              /* return; */
	    }
	}
      else
	{
	  tree->rates->t_prior_min[d->num] = tree->rates->t_prior_min[a->num];
	}
            
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      TIMES_Set_All_Node_Priors_Top_Down_S(d,d->v[i], result, tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl LOG_g_i(phydbl lmbd, phydbl t_slice_max, phydbl t_slice_min, phydbl t_prior_max, phydbl t_prior_min)
{
  phydbl result = 0.0;

  printf("\n. lmbd * t_prior_min = [%f] \n", lmbd * t_prior_min);
  if(lmbd * t_prior_min < -10.0)
    {
      phydbl K;
      K = -10.0 - lmbd * t_prior_min;
      printf("\n. K = [%f] \n", K);
      if(lmbd * t_prior_max + K > 700.0)
        {
          PhyML_Printf("\n. Please scale your calibration intervals. \n");
          PhyML_Printf("\n. Cannot calculate EXP(-lmbd * bound). \n");
          PhyML_Printf("\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
          Warn_And_Exit("\n");
        }
      else
        {
          printf("\n. slice_max [%f] slice_min [%f] t_prior_max [%f] t_prior_min [%f]\n", lmbd * t_slice_max + K, lmbd * t_slice_min + K, lmbd * t_prior_max + K, lmbd * t_prior_min + K);
          result = LOG(EXP(lmbd * t_slice_max + K) - EXP(lmbd * t_slice_min + K)) - LOG(EXP(lmbd * t_prior_max + K) - EXP(lmbd * t_prior_min + K));
        }
    }
  else 
    {
      result = LOG(EXP(lmbd * t_slice_max) - EXP(lmbd * t_slice_min)) - LOG(EXP(lmbd * t_prior_max) - EXP(lmbd * t_prior_min));
    }
  return(result);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
