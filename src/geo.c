/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


#include "geo.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int GEO_Main(int argc, char **argv)
{
  t_geo *t;
  int n_tax;
  t_tree *tree;
  int seed;
  int n_vars;
  phydbl *res;
  FILE *fp;
  int pid;
  char *s;

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  
  strcpy(s,"geo.out");
  pid = getpid();
  sprintf(s+strlen(s),".%d",pid);

  fp = fopen(s,"w");

  seed = time(NULL);
  /* seed = 1359056628; */
  printf("\n. Seed = %d",seed);
  srand(seed);

  t = GEO_Make_Geo_Basic();
  GEO_Init_Geo_Struct(t);

  /* t->tau        = 0.5; */
  /* t->lbda       = 0.2; */
  /* t->sigma      = 3.0; */

  t->tau        = Uni()*(5.0-0.1)   + 0.1;
  t->lbda       = Uni()*(5.0-0.1)   + 0.1;
  t->sigma      = Uni()*(10. - 0.1) + 0.1;

  t->ldscape_sz = 40;
  t->n_dim      = 2;
  n_tax         = 100;

  /* PhyML_Fprintf(fp,"\nSigmaTrue\t LbdaTrue\t TauTrue\tSigma5\t Sigma50\t Sigma95\t Lbda5\t Lbda50\t Lbda95\t Tau5\t Tau50\t Tau95\t"); */
  /* PhyML_Fprintf(fp,"\n"); */
  PhyML_Fprintf(fp,"%f\t %f\t %f\t",
                t->sigma,
                t->lbda,
                t->tau);


  GEO_Make_Geo_Complete(t->ldscape_sz,t->n_dim,n_tax,t);

  t->cov[0*t->n_dim+0] = t->sigma;
  t->cov[1*t->n_dim+1] = t->sigma;
  t->cov[0*t->n_dim+1] = 0.0;
  t->cov[1*t->n_dim+0] = 0.0;

  GEO_Simulate_Coordinates(t->ldscape_sz,t);

  tree = GEO_Simulate(t,n_tax);

  GEO_Get_Locations_Beneath(t,tree);

  GEO_Randomize_Locations(tree->n_root, 
                          tree->geo,
                          tree);
  GEO_Update_Occup(t,tree);
  GEO_Lk(t,tree);
  PhyML_Printf("\n. Init loglk: %f",tree->geo->c_lnL);

  /* Exit("\n"); */
  /* int c = 0; */
  /* do */
  /*   { */
  /*     GEO_Optimize_Sigma(t,tree); */
  /*     GEO_Optimize_Lambda(t,tree); */
  /*     GEO_Optimize_Tau(t,tree); */
  /*     printf("\n. Sigma = %f Lambda = %f Tau = %f Loglk = %f",tree->geo->sigma,tree->geo->lbda,tree->geo->tau,GEO_Lk(t,tree)); */
  /*     c++; */
  /*   } */
  /* while(c < 30); */


  tree->mcmc = MCMC_Make_MCMC_Struct();
  MCMC_Complete_MCMC(tree->mcmc,tree);

  tree->mcmc->io               = NULL;
  tree->mcmc->is               = NO;
  tree->mcmc->use_data         = YES;
  tree->mcmc->run              = 0;
  tree->mcmc->sample_interval  = 1E+3;
  tree->mcmc->chain_len        = 1E+6;
  tree->mcmc->chain_len_burnin = 1E+5;
  tree->mcmc->randomize        = YES;
  tree->mcmc->norm_freq        = 1E+3;
  tree->mcmc->max_tune         = 1.E+20;
  tree->mcmc->min_tune         = 1.E-10;
  tree->mcmc->print_every      = 2;
  tree->mcmc->is_burnin        = NO;
  tree->mcmc->nd_t_digits      = 1;


  t->tau   = 1.0;
  t->lbda  = 1.;
  t->sigma = 1.0;

  n_vars = 4;
  tree->mcmc->chain_len = 1.E+4;
  tree->mcmc->sample_interval = 50;


  GEO_Update_Occup(t,tree);
  GEO_Lk(t,tree);
  PhyML_Printf("\n. Init loglk: %f",tree->geo->c_lnL);


  res = (phydbl *)mCalloc(tree->mcmc->chain_len / tree->mcmc->sample_interval * n_vars,sizeof(phydbl));

  tree->mcmc->run = 0;
  do
    {
      MCMC_Geo_Lbda(tree);
      MCMC_Geo_Sigma(tree);
      MCMC_Geo_Tau(tree);
      MCMC_Geo_Loc(tree);
      
      if(tree->mcmc->run%tree->mcmc->sample_interval == 0)
        {
          PhyML_Printf("\n. Run %6d Sigma: %12f Lambda: %12f Tau: %12f LogLk: %12f",
                       tree->mcmc->run,
                       tree->geo->sigma,
                       tree->geo->lbda,
                       tree->geo->tau,
                       tree->geo->c_lnL);

          res[0 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->sigma; 
          res[1 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->lbda; 
          res[2 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->tau; 
          res[3 * tree->mcmc->chain_len / tree->mcmc->sample_interval +  tree->mcmc->run / tree->mcmc->sample_interval] = tree->geo->c_lnL; 
        }

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);
  
  
  PhyML_Fprintf(fp,"%f\t %f\t %f\t%f\t %f\t %f\t%f\t %f\t %f\n",
               Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.05),
               Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
               Quantile(res+0*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.95),

               Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.05),
               Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
               Quantile(res+1*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.95),

               Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.05),
               Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.50),
               Quantile(res+2*tree->mcmc->chain_len / tree->mcmc->sample_interval,tree->mcmc->run / tree->mcmc->sample_interval,0.95));

  Free(s);

  fclose(fp);

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_geo *GEO_Make_Geo_Basic()
{
  t_geo *t;
  t = (t_geo *)mCalloc(1,sizeof(t_geo));
  return(t);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Make_Geo_Complete(int ldscape_sz, int n_dim, int n_tax, t_geo *t)
{
  // F matrix
  t->f_mat = (phydbl *)mCalloc(ldscape_sz*ldscape_sz,sizeof(phydbl));

  // R matrix
  t->r_mat = (phydbl *)mCalloc(ldscape_sz*ldscape_sz,sizeof(phydbl));

  // Occupation vectors: one vector for each node
  t->occup = (int *)mCalloc((int)(2*n_tax-1)*ldscape_sz,sizeof(int)); 

  // Locations
  t->ldscape = (phydbl *)mCalloc((int)(ldscape_sz*n_dim),sizeof(phydbl));

  // Lineage locations
  t->loc = (int *)mCalloc((int)(2*n_tax-1),sizeof(int));

  // Sorted node heights
  t->sorted_nd = (t_node **)mCalloc((int)(2*n_tax-1),sizeof(t_node *));

  // Covariance matrix
  t->cov = (phydbl *)mCalloc((int)(n_dim*n_dim),sizeof(phydbl));

  // gives the location occupied beneath each node in the tree
  t->loc_beneath = (int *)mCalloc((int)(2*n_tax-1)*ldscape_sz,sizeof(int));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Geo(t_geo *t)
{
  Free(t->f_mat);
  Free(t->r_mat);
  Free(t->occup);
  Free(t->loc);
  Free(t->ldscape);
  Free(t->sorted_nd);
  Free(t->cov);
  Free(t->loc_beneath);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Update F matrix. Assume a diagonal covariance matrix.
void GEO_Update_Fmat(t_geo *t)
{
  phydbl *loc1, *loc2;
  int i,j,k;
  int err;
  phydbl lognloc;
  
  For(i,t->n_dim) t->cov[i*t->n_dim+i] = t->sigma; // Diagonal covariance matrix. Same variance in every direction

  lognloc = LOG(t->ldscape_sz);

  // Fill in F matrix;
  for(i=0;i<t->ldscape_sz;i++)
    {
      loc1 = t->ldscape + i*t->n_dim;

      for(j=i;j<t->ldscape_sz;j++)
        {
          loc2 = t->ldscape + j*t->n_dim;

          t->f_mat[i*t->ldscape_sz+j] = .0;

          // Calculate log(f(l_i,l_j)) - log(f(l_i,l_i) - log(m) 
          For(k,t->n_dim) t->f_mat[i*t->ldscape_sz+j] += Log_Dnorm(loc2[k],loc1[k],SQRT(t->cov[k*t->n_dim+k]),&err);
          t->f_mat[i*t->ldscape_sz+j] -= lognloc;
          For(k,t->n_dim) t->f_mat[i*t->ldscape_sz+j] -= Log_Dnorm(loc1[k],loc1[k],SQRT(t->cov[k*t->n_dim+k]),&err);

          // Take the exponential
          t->f_mat[i*t->ldscape_sz+j] = EXP(t->f_mat[i*t->ldscape_sz+j]);
          
          /* printf("\n. i=%d j=%d f=%f %f %f %f %f",i,j,t->f_mat[i*t->ldscape_sz+j],loc1[0],loc2[0],loc1[1],loc2[1]); */

          // Matrix is symmetric
          t->f_mat[j*t->ldscape_sz+i] = t->f_mat[i*t->ldscape_sz+j];
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Sort node heights from oldest to youngest age.
void GEO_Update_Sorted_Nd(t_geo *t, t_tree *tree)
{
  int i;
  int swap;
  t_node *buff;

  buff = NULL;

  For(i,2*tree->n_otu-1) t->sorted_nd[i] = tree->a_nodes[i];

  // Bubble sort of the node heights
  do
    {
      swap = NO;
      For(i,2*tree->n_otu-2) 
        {
          if(tree->rates->nd_t[t->sorted_nd[i+1]->num] < tree->rates->nd_t[t->sorted_nd[i]->num])
            {
              buff              = t->sorted_nd[i];
              t->sorted_nd[i]   = t->sorted_nd[i+1];
              t->sorted_nd[i+1] = buff;

              swap = YES;
            }
        }
    }
  while(swap == YES);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Update the set of vectors of occupation along the tree
void GEO_Update_Occup(t_geo *t, t_tree *tree)
{
  int i,j;

  GEO_Update_Sorted_Nd(t,tree);

  For(i,t->ldscape_sz*(2*tree->n_otu-1)) t->occup[i] = 0;

  t->occup[tree->n_root->num*t->ldscape_sz + t->loc[tree->n_root->num]]++;
  
  for(i=1;i<2*tree->n_otu-1;i++)
    {
      For(j,t->ldscape_sz) 
        {
          t->occup[t->sorted_nd[i]->num*t->ldscape_sz + j] = 
            t->occup[t->sorted_nd[i-1]->num*t->ldscape_sz + j];
        }

      if(t->sorted_nd[i-1]->tax == NO)
        t->occup[t->sorted_nd[i]->num*t->ldscape_sz + t->loc[t->sorted_nd[i]->num]]++;

      /* printf("\n ** %3d  ",t->sorted_nd[i]->num); */
      /* For(j,t->ldscape_sz) printf(" %d",t->occup[t->sorted_nd[i]->num*t->ldscape_sz + j]); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Calculate R mat at node n
void GEO_Update_Rmat(t_node *n, t_geo *t, t_tree *tree)
{
  int i,j;
  phydbl lbda_j;

  GEO_Update_Fmat(t);

  For(i,t->ldscape_sz)
    {
      For(j,t->ldscape_sz)
        {
          lbda_j = ((t->occup[n->num*t->ldscape_sz + j]>0) ? (t->lbda) : (1.0));
          t->r_mat[i*t->ldscape_sz+j] = t->f_mat[i*t->ldscape_sz+j] * lbda_j * t->tau;          
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl GEO_Lk(t_geo *t, t_tree *tree)
{
  int i;
  phydbl loglk;
  phydbl R;
  int dep,arr; // departure and arrival location indices;
  t_node *curr_n,*prev_n;
  phydbl sum;

  GEO_Update_Sorted_Nd(t,tree);
  GEO_Update_Occup(t,tree);     // Same here.

  prev_n = NULL;
  curr_n = NULL;
  loglk = .0;
  for(i=1;i<tree->n_otu-1;i++) // Consider all the time slices, from oldest to youngest. 
                             // Start at first node below root
    {
      prev_n = t->sorted_nd[i-1]; // node just above
      curr_n = t->sorted_nd[i];   // current node
      
      GEO_Update_Rmat(curr_n,t,tree); // NOTE: don't need to do that every time. Add check later.

      R = GEO_Total_Migration_Rate(curr_n,t); // Total migration rate calculated at node n

      dep = t->loc[curr_n->num]; // departure location
      arr =                      // arrival location
        (t->loc[curr_n->v[1]->num] == t->loc[curr_n->num] ? 
         t->loc[curr_n->v[2]->num] : 
         t->loc[curr_n->v[1]->num]);
      
      /* printf("\n%f\t%f", */
      /*        t->ldscape[arr*t->n_dim+0]-t->ldscape[dep*t->n_dim+0], */
      /*        t->ldscape[arr*t->n_dim+1]-t->ldscape[dep*t->n_dim+1]); */

      loglk -= R * FABS(tree->rates->nd_t[curr_n->num] - tree->rates->nd_t[prev_n->num]);
      loglk += LOG(t->r_mat[dep * t->ldscape_sz + arr]);

      /* printf("\n. R = %f r_mat = %f dt = %f loglk = %f", */
      /*        R, */
      /*        t->r_mat[dep * t->ldscape_sz + arr], */
      /*        FABS(t->sorted_nd[i] - t->sorted_nd[i-1]),loglk); */

    }


  // Likelihood for the first 'slice' (i.e., the part just below the root down to
  // the next node)
  GEO_Update_Rmat(tree->n_root,t,tree);

  loglk -= LOG(t->ldscape_sz); 
  dep = t->loc[tree->n_root->num];
  arr = 
    (t->loc[tree->n_root->num] != t->loc[tree->n_root->v[1]->num] ? 
     t->loc[tree->n_root->v[1]->num] :
     t->loc[tree->n_root->v[2]->num]);
  
  /* printf("\n %f %f",t->ldscape[dep],t->ldscape[arr]); */

  loglk += LOG(t->r_mat[dep * t->ldscape_sz + arr]);
    
  sum = .0;
  For(i,t->ldscape_sz) sum += t->r_mat[dep * t->ldscape_sz + i];
  
  loglk -= LOG(sum);

  tree->geo->c_lnL = loglk;

  return loglk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Init_Tloc_Tips(t_geo *t, t_tree *tree)
{
  int i;

  // TO DO
  For(i,tree->n_otu)
    {
      t->loc[i] = i;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Do not forget to call GEO_Update_Rmat (with node n) before calling this function
phydbl GEO_Total_Migration_Rate(t_node *n, t_geo *t)
{
  phydbl R;
  int i,j;

  R = .0;

  For(i,t->ldscape_sz)
    {
      For(j,t->ldscape_sz)
        {
          R += 
            t->r_mat[i*t->ldscape_sz + j] * 
            t->occup[n->num * t->ldscape_sz + i];
        }
    }

  return R;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Find the arrival location for the migration leaving from n
int GEO_Get_Arrival_Location(t_node *n, t_geo *t, t_tree *tree)
{
  int i;
  t_node *v1, *v2; // the two daughters of n

  v1 = v2 = NULL;

  For(i,3)
    {
      if(n->v[i] && n->v[i] != n->anc)
        {
          if(!v1) v1 = n->v[i];
          else    v2 = n->v[i];
        }
    }

  if(t->loc[v1->num] == t->loc[v2->num]) // Migrated to the same location as that of n
    {
      if(t->loc[n->num] != t->loc[v1->num])
        {
          PhyML_Printf("\n== Error detected in location labeling.");
          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");    
        }
      else
        return t->loc[n->num];
    }
  else // Migrated to a different spot
    {
      if((t->loc[v1->num] != t->loc[n->num]) && (t->loc[v2->num] != t->loc[n->num]))
        {
          PhyML_Printf("\n== Error detected in location labeling.");
          PhyML_Printf("\n== Err in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");    
        }
      else
        {
          if(t->loc[v1->num] == t->loc[n->num]) return t->loc[v2->num]; // v2 gets the new location, v1 inheritates from n
          else                                  return t->loc[v1->num]; // v1 gets the new location, v2 inheritates from n
        }
    }
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *GEO_Simulate(t_geo *t, int n_otu)
{
  t_tree *tree;
  int n_branching_nodes;
  t_node **branching_nodes; // vector of nodes available for branching out
  phydbl *p_branch; // p_branch[i]: proba that the i-th node in branching_nodes branches out (p_x vector in the article)
  phydbl *p_mig; // p_branch[i]: proba of migrating to location i from the location of the edge branching out (q_i vector in the article)
  int hit;
  phydbl time;
  int dep, arr;
  int i,j;
  phydbl sum;
  phydbl R;
  int *occup; // occupation vector. Updated as we move down the tree
  int nd_idx;
  t_node *buff_nd;
  phydbl buff_t;
  int buff_l;
  int swap;


  tree = Make_Tree_From_Scratch(n_otu,NULL);
  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,NULL,tree->n_otu);
  tree->n_root = tree->a_nodes[2*tree->n_otu-2]; // Set the root node to the last element in the list of nodes
  tree->geo = t;

  For(i,2*tree->n_otu-2) tree->rates->nd_t[i] = -1.;

  occup = (int *)mCalloc(t->ldscape_sz,sizeof(int));
  
  GEO_Update_Fmat(t);

  branching_nodes = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
  branching_nodes[0] = tree->n_root;
  n_branching_nodes  = 1;
  nd_idx = 0;

  
  p_branch = (phydbl *)mCalloc(tree->n_otu,sizeof(phydbl ));
  p_branch[0] = 1.0;

  p_mig = (phydbl *)mCalloc(t->ldscape_sz,sizeof(phydbl ));

  time = 0.0; // Time at the root node (this will be changed afterward)

  // Sample a location uniformly for the root
  t->loc[tree->n_root->num] = Rand_Int(0,t->ldscape_sz-1);
  
  // Update the occupancy vector
  occup[t->loc[tree->n_root->num]] = 1;

  dep = arr = -1;

 // total migration rate
  R = 0.0;
  For(i,t->ldscape_sz) 
    {
      R += 
        t->f_mat[t->loc[tree->n_root->num]*t->ldscape_sz+i] * 
        ((occup[i]>0) ? (t->lbda) : (1.0)) * 
        t->tau;
    }

  do
    {      
      // Select the node that branches out
      hit = Sample_i_With_Proba_pi(p_branch,n_branching_nodes);
      
      /* printf("\n. [%d] Select node %d (location %d)",n_branching_nodes,branching_nodes[hit]->num,t->loc[branching_nodes[hit]->num]); */

      // Set the time for the branching node
      tree->rates->nd_t[branching_nodes[hit]->num] = time;


      /* printf("\n. Set its time to %f",time); */

      // Select the destination location
      dep = t->loc[branching_nodes[hit]->num]; // Departure point
           
      sum = .0;
      For(i,t->ldscape_sz) // Total rate of migration out of departure point
        {
          p_mig[i] = 
            t->f_mat[dep*t->ldscape_sz+i] * 
            ((occup[i]>0) ? (t->lbda) : (1.0)) * 
            t->tau;

          sum += p_mig[i];
        }      
      For(i,t->ldscape_sz) p_mig[i] /= sum;

      arr = Sample_i_With_Proba_pi(p_mig,t->ldscape_sz);

      /* printf("\n. Migrate from %d [%5.2f,%5.2f] to %d [%5.2f,%5.2f]", */
      /*        dep, */
      /*        t->ldscape[dep*t->n_dim+0], */
      /*        t->ldscape[dep*t->n_dim+1], */
      /*        arr, */
      /*        t->ldscape[arr*t->n_dim+0], */
      /*        t->ldscape[arr*t->n_dim+1]); */

      /* printf("\n%f\t%f", */
      /*        t->ldscape[arr*t->n_dim+0]-t->ldscape[dep*t->n_dim+0], */
      /*        t->ldscape[arr*t->n_dim+1]-t->ldscape[dep*t->n_dim+1]); */
             


      // Update vector of occupation
      occup[arr]++;
      
      /* printf("\n. Remove %d. Add %d and %d",branching_nodes[hit]->num,tree->a_nodes[nd_idx]->num,tree->a_nodes[nd_idx+1]->num); */
      // Connect two new nodes to the node undergoing a branching event
      tree->a_nodes[nd_idx]->v[0]   = branching_nodes[hit];
      tree->a_nodes[nd_idx+1]->v[0] = branching_nodes[hit];
      branching_nodes[hit]->v[1] = tree->a_nodes[nd_idx];
      branching_nodes[hit]->v[2] = tree->a_nodes[nd_idx+1];

      // update branching_nodes vector. Element 'hit' is being replaced so that the corresponding node can no longer branch out
      branching_nodes[hit]                 = tree->a_nodes[nd_idx];
      branching_nodes[n_branching_nodes]   = tree->a_nodes[nd_idx+1];

      // Update t_loc vector.
      t->loc[tree->a_nodes[nd_idx]->num]   = dep;
      t->loc[tree->a_nodes[nd_idx+1]->num] = arr;

      // Update total migration rate 
      R = .0;
      For(i,t->ldscape_sz)
        {
          if(occup[i] > 0)
            {
              For(j,t->ldscape_sz)
                {
                  R += 
                    occup[i] *
                    t->f_mat[i*t->ldscape_sz+j] * 
                    ((occup[j]>0) ? (t->lbda) : (1.0)) *
                    t->tau;
                }
            }
        }

      // Set the time until next branching event
      time = time + Rexp(R);
    
      // Update p_branch vector
      For(i,n_branching_nodes+1)
        {
          dep = t->loc[branching_nodes[i]->num];
          p_branch[i] = 0.0;
          For(j,t->ldscape_sz)
            {
              p_branch[i] +=
                t->f_mat[dep*t->ldscape_sz+j] * 
                ((occup[j]>0) ? (t->lbda) : (1.0)) * 
                t->tau / R;

              /* printf("\n. %f %f %f %f", */
              /*        R, */
              /*        t->f_mat[dep*t->ldscape_sz+j], */
              /*        ((occup[j]>0) ? (t->lbda) : (1.0)), */
              /*        t->tau); */
            }
          /* printf("\n. %f ",p_branch[i]); */
        }

              
      // Increase the number of branching nodes by one (you just added 2 new and removed 1 old)
      n_branching_nodes++;
      nd_idx += 2;

    }
  while(n_branching_nodes < n_otu);

  printf("\n");
  For(i,t->ldscape_sz)
    printf("%d,",occup[i]);
  printf("\n");


  // Set the times at the tips
  For(i,2*tree->n_otu-1) if(tree->rates->nd_t[i] < 0.0) tree->rates->nd_t[i] = time;
  
  // Reverse time scale
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] -= time;
  /* For(i,2*tree->n_otu-1) tree->rates->nd_t[i] = FABS(tree->rates->nd_t[i]); */

  //  Bubble sort to put all the tips at the top of the tree->a_nodes array
  do
    {
      swap = NO;
      For(i,2*tree->n_otu-2)
        {
          if(!tree->a_nodes[i+1]->v[1] && tree->a_nodes[i]->v[1])
            {
              buff_nd            = tree->a_nodes[i+1];
              tree->a_nodes[i+1] = tree->a_nodes[i];
              tree->a_nodes[i]   = buff_nd;

              buff_t                 = tree->rates->nd_t[i+1];
              tree->rates->nd_t[i+1] = tree->rates->nd_t[i];
              tree->rates->nd_t[i]   = buff_t;

              buff_l      = t->loc[i+1];
              t->loc[i+1] = t->loc[i];
              t->loc[i]   = buff_l;

              swap = YES;
            }
        }
    }
  while(swap == YES);


  // The rest below is just bookkeeping...


  For(i,2*tree->n_otu-1) tree->a_nodes[i]->num = i;
  For(i,2*tree->n_otu-1) 
    {
      if(i < tree->n_otu) tree->a_nodes[i]->tax = YES;
      else                tree->a_nodes[i]->tax = NO;
    }

  /* printf("\n++++++++++++++++++\n"); */
  /* For(i,2*tree->n_otu-1) */
  /*   { */
  /*     printf("\n. Node %3d [%p] anc:%3d v1:%3d v2:%3d time: %f", */
  /*            tree->a_nodes[i]->num, */
  /*            (void *)tree->a_nodes[i], */
  /*            tree->a_nodes[i]->v[0] ? tree->a_nodes[i]->v[0]->num : -1, */
  /*            tree->a_nodes[i]->v[1] ? tree->a_nodes[i]->v[1]->num : -1, */
  /*            tree->a_nodes[i]->v[2] ? tree->a_nodes[i]->v[2]->num : -1, */
  /*            tree->rates->nd_t[i]);              */
    /* } */


  For(i,tree->n_otu) 
    {
      if(!tree->a_nodes[i]->name) tree->a_nodes[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(tree->a_nodes[i]->name,"x");
      sprintf(tree->a_nodes[i]->name+1,"%d",i);      
    }

  tree->n_root->v[1]->v[0] = tree->n_root->v[2];
  tree->n_root->v[2]->v[0] = tree->n_root->v[1];
  
  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->a_nodes[0],tree->a_nodes[0]->v[0],tree);

  tree->e_root = tree->n_root->v[1]->b[0];

  For(i,2*tree->n_otu-3)
    {
      tree->a_edges[i]->l->v = FABS(tree->rates->nd_t[tree->a_edges[i]->left->num] - 
                                    tree->rates->nd_t[tree->a_edges[i]->rght->num]);
    }

  tree->e_root->l->v =
    FABS(tree->rates->nd_t[tree->n_root->v[1]->num] -
         tree->rates->nd_t[tree->n_root->num]) +
    FABS(tree->rates->nd_t[tree->n_root->v[2]->num] -
         tree->rates->nd_t[tree->n_root->num]);

  tree->n_root_pos = 
    FABS(tree->rates->nd_t[tree->n_root->v[2]->num] -
         tree->rates->nd_t[tree->n_root->num]) / tree->e_root->l->v;

  /* printf("\n. %s ",Write_Tree(tree,NO)); */

  DR_Draw_Tree("essai.ps",tree);

  /* For(i,tree->n_otu) */
  /*   printf("\n. %4s %4d [%5.2f %5.2f]",tree->a_nodes[i]->name, */
  /*          t->loc[i], */
  /*          t->ldscape[t->loc[i]*t->n_dim+0], */
  /*          t->ldscape[t->loc[i]*t->n_dim+1]); */

  
  Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);		

  Free(branching_nodes);
  Free(p_branch);
  Free(p_mig);

  return(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Simualte n coordinates (in 2D)
void GEO_Simulate_Coordinates(int n, t_geo *t)
{
  int i;
  phydbl width;

  width = 5.;

  For(i,n)
    {
      /* t->ldscape[i*t->n_dim+0] = Rnorm(0.0,3.); */
      /* t->ldscape[i*t->n_dim+1] = Rnorm(0.0,3.); */
      t->ldscape[i*t->n_dim+0] = -width/2. + Uni()*width;
      t->ldscape[i*t->n_dim+1] = -width/2. + Uni()*width;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Optimize_Sigma(t_geo *t, t_tree *tree)
{
  Generic_Brent_Lk(&(t->sigma),
                   t->min_sigma,
                   t->max_sigma,
                   1.E-5,
                   1000,
                   NO,
                   GEO_Wrap_Lk,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Optimize_Lambda(t_geo *t, t_tree *tree)
{
  Generic_Brent_Lk(&(t->lbda),
                   t->min_lbda,
                   t->max_lbda,
                   1.E-5,
                   1000,
                   NO,
                   GEO_Wrap_Lk,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Optimize_Tau(t_geo *t, t_tree *tree)
{
  Generic_Brent_Lk(&(t->tau),
                   t->min_tau,
                   t->max_tau,
                   1.E-5,
                   1000,
                   NO,
                   GEO_Wrap_Lk,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl GEO_Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return GEO_Lk(tree->geo,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Init_Geo_Struct(t_geo *t)
{
  t->c_lnL      = UNLIKELY;

  t->sigma      = 1.0;
  t->min_sigma  = 1.E-3;
  t->max_sigma  = 1.E+3;

  t->lbda       = 1.0;
  t->min_lbda   = 1.E-3;
  t->max_lbda   = 1.E+3;
  
  t->tau        = 1.0;
  t->min_tau    = 1.E-3;
  t->max_tau    = 1.E+3;

  t->tau        = 1.0;

  t->n_dim      = -1;
  t->ldscape_sz = -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Randomize_Locations_Pre(t_node *n, t_geo *t, t_tree *tree)
{

  if(n->tax == YES) return;
  else
    {
      int i;
      
      GEO_Randomize_Locations(n,t,tree);

      For(i,3)
        if(n->v[i] != n->anc && n->b[i] != tree->e_root) 
          GEO_Randomize_Locations_Pre(n->v[i],t,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Randomize_Locations(t_node *n, t_geo *t, t_tree *tree)
{
  t_node *v1, *v2;
  int i;
  phydbl *probs; // vector of probability of picking each location
  phydbl sum;
  
  probs = (phydbl *)mCalloc(t->ldscape_sz,sizeof(phydbl));

  v1 = v2 = NULL;
  For(i,3)
    {
      if(n->v[i] != n->anc && n->b[i] != tree->e_root)
        {
          if(!v1) v1 = n->v[i];
          else    v2 = n->v[i];
        }
    }

  if(t->loc[v1->num] != t->loc[n->num] && v1->tax == NO)
    {
      sum = 0.0;
      For(i,t->ldscape_sz) sum += t->loc_beneath[v1->num * t->ldscape_sz + i];
      For(i,t->ldscape_sz) probs[i] = t->loc_beneath[v1->num * t->ldscape_sz + i]/sum;
      
      t->loc[v1->num] = Sample_i_With_Proba_pi(probs,t->ldscape_sz);      
    }
  else if(t->loc[v2->num] != t->loc[n->num] && v2->tax == NO)
    {
      sum = 0.0;
      For(i,t->ldscape_sz) sum += t->loc_beneath[v2->num * t->ldscape_sz + i];
      For(i,t->ldscape_sz) probs[i] = t->loc_beneath[v2->num * t->ldscape_sz + i]/sum;
      
      Sample_i_With_Proba_pi(probs,t->ldscape_sz);

      t->loc[v2->num] = Sample_i_With_Proba_pi(probs,t->ldscape_sz);      
    }

  Free(probs);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Get_Locations_Beneath(t_geo *t, t_tree *tree)
{
  int i;

  GEO_Get_Locations_Beneath_Post(tree->n_root,tree->n_root->v[1],t,tree);
  GEO_Get_Locations_Beneath_Post(tree->n_root,tree->n_root->v[2],t,tree);

  For(i,t->ldscape_sz)
    t->loc_beneath[tree->n_root->num*t->ldscape_sz+i] =
    t->loc_beneath[tree->n_root->v[1]->num*t->ldscape_sz+i] +
    t->loc_beneath[tree->n_root->v[2]->num*t->ldscape_sz+i] ;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void GEO_Get_Locations_Beneath_Post(t_node *a, t_node *d, t_geo *t, t_tree *tree)
{

  if(d->tax) 
    {
      t->loc_beneath[d->num*t->ldscape_sz+t->loc[d->num]] = 1;
      return;
    }
  else
    {
      int i;
      t_node *v1, *v2;

      For(i,3) 
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              GEO_Get_Locations_Beneath_Post(d,d->v[i],t,tree);
            }
        }
          

      v1 = v2 = NULL;
      For(i,3) 
        {
          if(d->v[i] != a && d->b[i] != tree->e_root)
            {
              if(!v1) v1 = d->v[i];
              else    v2 = d->v[i];
            }
        }
          
      For(i,t->ldscape_sz) 
        t->loc_beneath[ d->num*t->ldscape_sz+i] = 
        t->loc_beneath[v1->num*t->ldscape_sz+i] + 
        t->loc_beneath[v2->num*t->ldscape_sz+i] ;
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


