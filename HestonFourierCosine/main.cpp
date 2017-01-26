#include <iostream>
#include <stdio.h>
#include <math.h>
#include "de.h"
#include "dtime.h"
#include <ctime>
#include <nlopt.h>
#include "wrapper.h"

double fa_minbound[5] = {1e-8, 1e-8, 1e-8, -1.0+1e-8, 1e-8};
double fa_maxbound[5] = {5.0-1e-8, 1.0-1e-8, 1.0-1e-8, 1.0-1e-8, 1.0-1e-8};

double dtime();
t_pop devol(FILE *, FILE *, double[], double[]);

int main()
{

     	FILE *Fp_in;
     	FILE *Fp_out;
	double tcos = 0.0;
	t_pop best;
	nlopt_opt opt;
	double minf=0.0, x[5];// = {4.7244004508, 0.3343881602, 0.8395414554, -0.7663721088, 0.7069578399};

	Fp_in   = fopen("in.dat","r");
     	Fp_out  = fopen("out.dat","w");  
	 
	if(taking_input())
	{

		 double t1 = dtime();
		 //---optimization--------------------------------------
		 best = devol( Fp_in, Fp_out, fa_minbound, fa_maxbound );
		 double t2 = dtime();
		 
		 tcos += t2-t1;
		 
		 std::cout<<"DEoptim Time = "<<tcos<<" sec"<<std::endl;
		 
		 fprintf(Fp_out, "RMSE = %12.10f\n", best.fa_cost[0]);
		 fprintf(Fp_out,"\n******** best vector ********\n");
		 for (int i=0; i < MAXDIM-1; i++)
		 {
			 fprintf(Fp_out,"best_vector[%d]=%12.10f\n", i, best.fa_vector[i]);
			 x[i] = best.fa_vector[i];
		 }
		 fprintf(Fp_out,"\n");
		 
		 double t3 = dtime();
		 
		 opt = nlopt_create(NLOPT_LN_COBYLA,5);
		 nlopt_set_lower_bounds(opt, fa_minbound);
		 nlopt_set_upper_bounds(opt, fa_maxbound);
		 nlopt_set_min_objective(opt, c_main, NULL);
		 nlopt_set_xtol_rel(opt, 1e-4);

		 nlopt_set_maxeval(opt,100);
		 
		 if (nlopt_optimize(opt, x, &minf) < 0) 
		 {
			printf("nlopt failed!\n");
		 }
		 else 
		 {
			//printf("found minimum at  f(%g,%g) = %0.10g\n", x[0], x[1], minf);
		 }
		 
		 double t4 = dtime(); 

		 tcos +=  t4-t3;
		 std::cout<<"Nlopt Time = "<<t4-t3<<" sec"<<std::endl;
		 
		 fprintf(Fp_out, "NLOPT RMSE = %f\n", minf);
		 fprintf(Fp_out,"\n******** best vector ********\n");
		 
		 for (int i=0; i < MAXDIM-1; i++)
		 {
			 fprintf(Fp_out,"best_vector[%d]=%12.10f\n", i, x[i]);
		 }
		 
		 std::cout<<"Total Time = "<<tcos<<" sec"<<std::endl;
		 free();
		 
		 fclose(Fp_in);
		 fclose(Fp_out);
	}
	else
	   	 printf("File Not Found\n");

	return 0;
 }
