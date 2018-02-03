/*
* Implements the Differential Evolution (DE) method
* for continuous function optimization by using MKL, Heston
* Fourier Cosine, compiler vectorisation with #pragma omp simd 
* directives, plus OpenMP parallelisation
*
* author: Matthew Dixon and Xin Zhang,  based on previous code written by 
*         Rainer Storn
*/

   
#define BOUND_CONSTR    //If defined the bounds fa_minbound[] and fa_maxbound[]
                        //are not only used for initializing the vector population
                        //but also used to keep the population within these bounds  
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "random.h"
#include "de.h"
#include <ctime>

extern double gfa_bound[21];

int   gi_gen;             // generation counter
int   gi_strategy;        // chooses DE-strategy
long  gl_nfeval;          // number of function evaluations     
int   gi_D;               // Dimension of parameter vector
int   gi_NP;              // Number of population members
int   gi_genmax;          // Maximum number of generations
t_pop gta_pop[2*MAXPOP];  // the two populations are put into one array side by side.      
t_pop gt_best;            // current best population member
t_pop *gpta_old, *gpta_new, *gpta_swap;


//---------------------General functions------------------------------

t_pop evaluate(int i_D, t_pop t_tmp, long *l_nfeval, t_pop *tpa_array, int i_NP);
int    left_vector_wins(t_pop trial, t_pop target);
void   devol(FILE *Fp_in, FILE *Fp_out);
void   sort (t_pop ary[], int len);
void   assigna2b(int D, double a[], double b[]);
void   permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid);


//-----------------------Main function--------------------------------                                                                             

int main()
 {

     FILE *Fp_in;
     FILE *Fp_out;
     double tcos = 0.0;
     Fp_in   = fopen("in.dat","r");
     Fp_out  = fopen("out.dat","w");             
     clock_t t1, t2;
     t1 = clock();
     //------------------optimization---------------------------------
	 devol(Fp_in,Fp_out);
	 t2 = clock();
         tcos += (((double)t2 - (double)t1) / CLOCKS_PER_SEC );
	 std::cout<<"Elapsed Time = "<<tcos<<" sec"<<std::endl;	 
	 fclose(Fp_in);
	 fclose(Fp_out);
	 return 0;
 }


//---------Assigns i_D-dimensional vector fa_a to vector f_b----------

/*
*  Parameters: i_D       size of vectors
*              fa_a[]    source vector
*              fa_b[]    destination vector
*/
void  assigna2b(int i_D, double fa_a[], double fa_b[])                                                                                                                                          {
   int j;
   for (j=0; j<i_D; j++)
   {
      fa_b[j] = fa_a[j];
   }
}

                                                                                                                                   
//----------Performs Differential Evolution optimization--------------                 
/*                                                                                                                                                                                 
*  Description:    gi_strategy  (I) 1 --> DE/rand/1:
*                                         the classical version of DE.             
*                                   2 --> DE/local-to-best/1:
*                                         a version which has been used by quite a number
*                                         of scientists. Attempts a balance between robustness
*                                         and fast convergence.
*                                   3 --> DE/best/1 with jitter:
*                                         taylored for small population sizes and fast convergence.
*                                         Dimensionality should not be too high.           
*                                   4 --> DE/rand/1 with per-vector-dither:
*                                         classical DE with dither to become even more robust.
*                                   5 --> DE/rand/1 with per-generation-dither:
*                                         classical DE with dither to become even more robust.
*                                         Choosing f_weight = 0.3 is a good start here.
*                                   6 --> DE/rand/1 either-or-algorithm:
*                                         Alternates between differential mutation and three-point-
*                                         recombination.           
* 
*                  gi_D         (I)    number of parameters 
*                  gl_nfeval   (I/O)   counter for function evaluations                                               
*                  gi_gen      (I/O)   number of generations  
*                  gt_best     (I/O)   best vector
*                                                 
*  Parameters:     Fp_in    (I)     pointer to input file
*                  Fp_out   (I)     pointer to output file   
*/                                                                  
void devol(FILE *Fp_in, FILE *Fp_out)
{
#define URN_DEPTH   5   //4 + one index to avoid

//--------------------Variable declarations---------------------------

   std::cout.precision(16);
   int   i, j, k;                 // counting variables                 
   int   i_r1, i_r2, i_r3, i_r4;  // placeholders for random indexes    
   int   i_refresh;               // refresh rate of screen output      
   int   i_genmax, i_seed, i_bs_flag;
   int   ia_urn2[URN_DEPTH];

   double fa_minbound[MAXDIM]={1e-8, 1e-8, 1e-8, -1.0 + 1e-8, 1e-8};
   double fa_maxbound[MAXDIM]={5 - 1e-8, 1 - 1e-8, 1 - 1e-8, 1 - 1e-8, 1- 1e-8};
   t_pop t_tmp, t_bestit;
#ifdef BOUND_CONSTR
   t_pop t_origin;
#endif//BOUND_CONSTR
   double f_weight, f_jitter, f_dither;
   double f_cross;


//--------------Initialization of annealing parameters----------------

 fscanf(Fp_in,"%d",&gi_strategy);      //---choice of strategy-----------------
 fscanf(Fp_in,"%d",&i_genmax);         //---maximum number of generations------
 gi_genmax = i_genmax;
 fscanf(Fp_in,"%d",&gi_D);             //---number of parameters---------------
 if (gi_D > MAXDIM)
 {
	 printf("Error! too many parameters\n");
	 return;
 }
 fscanf(Fp_in,"%d",&gi_NP);             //---population size.-------------------
 if (gi_NP > MAXPOP)
 {
	 printf("Error! too many points\n");
	 return;
 }
 fscanf(Fp_in,"%lf",&f_weight);         //---weight factor----------------------
 fscanf(Fp_in,"%lf",&f_cross);          //---crossing over factor---------------
 fscanf(Fp_in,"%d",&i_seed);            //---random seed------------------------
 fscanf(Fp_in,"%d",&i_bs_flag);         //---if TRUE: best of parent+child selection--------
                                        //---if FALSE: DE standard tournament selection----


//------------Initialize random number generator-------------------

 sgenrand((unsigned long)i_seed);// for Mersenne Twister
 gl_nfeval    =  0;              // reset number of function evaluations 

//--------------------Initialization-----------------------------
   for (j=0; j<gi_D; j++)
   {
	 gta_pop[0].fa_vector[j] = fa_minbound[j]+genrand()*(fa_maxbound[j] - fa_minbound[j]);
   }
   
   gta_pop[0]      = evaluate(gi_D,gta_pop[0],&gl_nfeval,&gta_pop[0],gi_NP);
   gt_best  = gta_pop[0];
      
   for (i=1; i<gi_NP; i++)
   {
      for (j=0; j<gi_D; j++)
      {
	     gta_pop[i].fa_vector[j] = fa_minbound[j]+genrand()*(fa_maxbound[j] - fa_minbound[j]);
      }
      gta_pop[i] = evaluate(gi_D,gta_pop[i],&gl_nfeval,&gta_pop[0],gi_NP);

      if(left_vector_wins(gta_pop[i],gt_best) == TRUE)
	{
          gt_best = gta_pop[i];
	}
   }

   t_bestit  = gt_best;

   //---assign pointers to current ("old") and new population---

   gpta_old = &gta_pop[0];
   gpta_new = &gta_pop[gi_NP];

  
//---------------------Iteration loop-------------------------
   gi_gen = 0;
   //Note that kbhit() needs conio.h which is not always available under Unix.
   while ((gi_gen < i_genmax))// && (kbhit() == 0))// && (gt_best.fa_cost[0] > VTR))
   {
          gi_gen++;
	  printf("Max Gen: %d\n",gi_gen);
	  //-------computer dithering factor (if needed)-------------
	  f_dither = f_weight + genrand()*(1.0 - f_weight);

      //--------start of loop through ensemble--------------
      for (i=0; i<gi_NP; i++)           
      {
		permute(ia_urn2,URN_DEPTH,gi_NP,i); //Pick 4 random and distinct
		i_r1 = ia_urn2[1];                 //population members
		i_r2 = ia_urn2[2];
		i_r3 = ia_urn2[3];
		i_r4 = ia_urn2[4];

        //------------Choice of strategy-----------------
		//--------classical strategy DE/rand/1/bin----------
		if (gi_strategy == 1)
		{
                 assigna2b(gi_D,gpta_old[i].fa_vector,t_tmp.fa_vector);
	         j = (int)(genrand()*gi_D); // random parameter         
	         k = 0;
	         do
		   {                            // add fluctuation to random target 
	             t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + f_weight*(gpta_old[i_r2].fa_vector[j]-gpta_old[i_r3].fa_vector[j]);
	             j = (j+1)%gi_D;
	             k++;
		   }while((genrand() < f_cross) && (k < gi_D));
#ifdef BOUND_CONSTR
		   assigna2b(gi_D,gpta_old[i_r1].fa_vector,t_origin.fa_vector);
#endif//BOUND_CONSTR
		}
		//-------------DE/local-to-best/1/bin----------------
	       else if (gi_strategy == 2)
		{
                  assigna2b(gi_D,gpta_old[i].fa_vector,t_tmp.fa_vector);
	          j = (int)(genrand()*gi_D);    // random parameter         
	          k = 0;
	          do
		   {                            // add fluctuation to random target 
	             t_tmp.fa_vector[j] = t_tmp.fa_vector[j] + f_weight*(t_bestit.fa_vector[j] - t_tmp.fa_vector[j]) + 
				                                        f_weight*(gpta_old[i_r2].fa_vector[j]-gpta_old[i_r3].fa_vector[j]);
	             j = (j+1)%gi_D;
	             k++;
		   }while((genrand() < f_cross) && (k < gi_D));
#ifdef BOUND_CONSTR
		   assigna2b(gi_D,t_tmp.fa_vector,t_origin.fa_vector);
#endif//BOUND_CONSTR
		}
		//-------------DE/best/1/bin with jitter----------------
	      else if (gi_strategy == 3)
		{
                  assigna2b(gi_D,gpta_old[i].fa_vector,t_tmp.fa_vector);
	          j = (int)(genrand()*gi_D); // random parameter         
	          k = 0;
	          do
		    {                          // add fluctuation to random target 
			  f_jitter = (0.0001*genrand()+f_weight);
	                  t_tmp.fa_vector[j] = t_bestit.fa_vector[j] + f_jitter*(gpta_old[i_r1].fa_vector[j]-gpta_old[i_r2].fa_vector[j]);
	                  j = (j+1)%gi_D;
	                  k++;
		   }while((genrand() < f_cross) && (k < gi_D));
#ifdef BOUND_CONSTR
		   assigna2b(gi_D,t_tmp.fa_vector,t_origin.fa_vector);
#endif//BOUND_CONSTR
		}
		//----------DE/rand/1/bin with per-vector-dither--------------
	      else if (gi_strategy == 4)
		{
                  assigna2b(gi_D,gpta_old[i].fa_vector,t_tmp.fa_vector);
	          j = (int)(genrand()*gi_D); // random parameter         
	          k = 0;
	          do
		    {                            // add fluctuation to random target 
	               t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + 
			                       (f_weight + genrand()*(1.0 - f_weight))*
								   (gpta_old[i_r2].fa_vector[j]-gpta_old[i_r3].fa_vector[j]);
	               j = (j+1)%gi_D;
	               k++;
		   }while((genrand() < f_cross) && (k < gi_D));
#ifdef BOUND_CONSTR
		   assigna2b(gi_D,t_tmp.fa_vector,t_origin.fa_vector);
#endif//BOUND_CONSTR
		}
		//------------DE/rand/1/bin with per-generation-dither-------------
	      else if (gi_strategy == 5)
               {
                 assigna2b(gi_D,gpta_old[i].fa_vector,t_tmp.fa_vector);
	         j = (int)(genrand()*gi_D); // random parameter         
	         k = 0;
	         do
		   {                            // add fluctuation to random target 
	          t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + f_dither*(gpta_old[i_r2].fa_vector[j]-gpta_old[i_r3].fa_vector[j]);
	          j = (j+1)%gi_D;
	          k++;
		   }while((genrand() < f_cross) && (k < gi_D));
#ifdef BOUND_CONSTR
		   assigna2b(gi_D,t_tmp.fa_vector,t_origin.fa_vector);
#endif//BOUND_CONSTR
		}
		//----------variation to DE/rand/1/bin: either-or-algorithm---------------
	    else
		{
                  assigna2b(gi_D,gpta_old[i].fa_vector,t_tmp.fa_vector);
	          j = (int)(genrand()*gi_D); // random parameter         
	          k = 0;
		   if (genrand() < 0.5)      //Pmu = 0.5
		   {                         //differential mutation
	             do
		       {                            // add fluctuation to random target 
	                 t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + f_weight*(gpta_old[i_r2].fa_vector[j]-gpta_old[i_r3].fa_vector[j]);
	                 j = (j+1)%gi_D;
	                 k++;
		       }while((genrand() < f_cross) && (k < gi_D));
		   }
		   else
		   {//recombination with K = 0.5*(F+1) --> F-K-Rule
	          do
	            {                            // add fluctuation to random target 
	              t_tmp.fa_vector[j] = gpta_old[i_r1].fa_vector[j] + 0.5*(f_weight+1.0)*
					                 (gpta_old[i_r2].fa_vector[j]+gpta_old[i_r3].fa_vector[j] -
									2*gpta_old[i_r1].fa_vector[j]);
	              j = (j+1)%gi_D;
	              k++;
		    }while((genrand() < f_cross) && (k < gi_D));
		   }
#ifdef BOUND_CONSTR
		   assigna2b(gi_D,gpta_old[i_r1].fa_vector,t_origin.fa_vector);
#endif//BOUND_CONSTR		
		}//end if (gi_strategy ...

#ifdef BOUND_CONSTR
      for (j=0; j<gi_D; j++) //----boundary constraints via random reinitialization-------
      {                      //--------------and bounce back--------------------
         if (t_tmp.fa_vector[j] < fa_minbound[j])
		 {
            t_tmp.fa_vector[j] = fa_minbound[j]+genrand()*(t_origin.fa_vector[j] - fa_minbound[j]);
		 }
         if (t_tmp.fa_vector[j] > fa_maxbound[j])
		 {
            t_tmp.fa_vector[j] = fa_maxbound[j]+genrand()*(t_origin.fa_vector[j] - fa_maxbound[j]);
		 }
      }
#endif//BOUND_CONSTR

        //---------Trial mutation now in t_tmp-----------------
	    t_tmp = evaluate(gi_D,t_tmp,&gl_nfeval,&gpta_old[0],gi_NP);  // Evaluate mutant in t_tmp[] 

if (i_bs_flag == TRUE)
{
        gpta_new[i]=t_tmp; //save new vector, selection will come later
}
else
{
	if (left_vector_wins(t_tmp,gpta_old[i]) == TRUE)
	   {
	        gpta_new[i]=t_tmp;                 // replace target with mutant 
		if (left_vector_wins(t_tmp,gt_best) == TRUE)// Was this a new minimum? 
		   {                               // if so...
		      gt_best = t_tmp;             // store best member so far  
		   }                               // If mutant fails the test...
	   }                                       // go to next the configuration 
	else
	   {
	      gpta_new[i]=gpta_old[i];             // replace target with old value 
	   }
}                                                  //if (i_bs_flag == TRUE)
     }                                            //end for (i=0; i<gi_NP; i++)
					    // End mutation loop through pop. 

if (i_bs_flag == TRUE)
{
      sort(gpta_old, 2*gi_NP); //sort array of parents + children
	  gt_best = gpta_old[0];
}
else
{
      gpta_swap = gpta_old;
      gpta_old  = gpta_new;
      gpta_new  = gpta_swap;
}//if (i_bs_flag == TRUE)
	  t_bestit = gt_best;
	  
	  printf("%ld   %lf   %lf\n",gl_nfeval, gt_best.fa_cost[0], gt_best.fa_constraint[0]);

   }//end while ((gi_gen < i_genmax) && (gf_emin > MINI))
   
   fprintf(Fp_out,"RMSE = %12.10f\n", gt_best.fa_cost[0]);
   fprintf(Fp_out,"\n******** best vector ********\n", i, gt_best.fa_vector[i]);
   for (i=0; i<gi_D; i++)
   {
	   fprintf(Fp_out,"best_vector[%d]=%12.10f\n", i, gt_best.fa_vector[i]);
   }
   
}


//----------------Generate i_urn2_depth random indices-------------------
/*                                                                                                                                
* Description:    Generates i_urn2_depth random indices ex [0, i_NP-1]
*                 which are all distinct. This is done by using a 
*                 permutation algorithm called the "urn algorithm"
*                 which goes back to C.L.Robinson.                                                               
* Parameters:     ia_urn2       (O)    array containing the random indices
*                 i_urn2_depth  (I)    number of random indices (avoided index included)
*                 i_NP          (I)    range of indices is [0, i_NP-1]
*                 i_avoid       (I)    is the index to avoid and is located in
*                                      ia_urn2[0].   
*                                                                  
* Preconditions:  # Make sure that ia_urn2[] has a length of i_urn2_depth.
*                 # i_urn2_depth must be smaller than i_NP.                     
*                                                                  
* Postconditions: # the index to be avoided is in ia_urn2[0], so fetch the
*                   indices from ia_urn2[i], i = 1, 2, 3, ..., i_urn2_depth.                                                                                                               
*/

void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid)
{
	int  i, k, i_urn1, i_urn2;
	int  ia_urn1[MAXPOP] = {0};      //urn holding all indices
	k      = i_NP;
	i_urn1 = 0; 
	i_urn2 = 0;
	for (i=0; i<i_NP; i++) ia_urn1[i] = i; //initialize urn1
	i_urn1 = i_avoid;                  //get rid of the index to be avoided and place it in position 0.
	while (k >= i_NP-i_urn2_depth)     //i_urn2_depth is the amount of indices wanted (must be <= NP) 
	{
	   ia_urn2[i_urn2] = ia_urn1[i_urn1];      //move it into urn2
	   ia_urn1[i_urn1] = ia_urn1[k-1]; //move highest index to fill gap
	   k = k-1;                        //reduce number of accessible indices
	   i_urn2 = i_urn2 + 1;            //next position in urn2
           i_urn1 = (int)(genrand()*k);    //choose a random index 
        }
}
