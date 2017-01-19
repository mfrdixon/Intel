/* Sample code implementing the cost function evaluation
*  by using selected DE criterion.
*
*  author: Matthew Dixon and Xin Zhang
*  Date: Nov.30 2016
*/

#include "de.h"
#include<stdio.h>
#include<stdlib.h>
#include "wrapper.h"
#include <iostream>

//------global variable for the objective function--------------
float gfa_bound[21] ={0,1.2,1.88,3.119,6.06879,
		  11.25312,20.93868,38.99973,
		  72.66066687999998,135.385869312,
		  252.2654194687999, 470.0511374131198,
		  875.857310322687,  1632.00640736133,
		  3040.958067344504, 5666.29295426548,
		  10558.14502289265, 19673.25510067688,
		  36657.66721873185, 68305.14622427958,
		  127274.6837195391};
		  
		 
//----------------------objective function------------------------------
/*
* Description:   Evaluates the actual cost function (objective function)
*                which in this case evaluates the Chebychev fitting problem
*
* Parameters:    i_D        number of parameters
*                t_tmp      parameter vector
*                l_nfeval   counter for function evaluations
*                tpa_array  pointer to current population (not needed here)
*                i_NP       number of population members (not needed here)
* 
* Return value:  TRUE if trial vector wins, FALSE otherwise
*
*/
t_pop evaluate(int i_D, t_pop t_tmp, long *l_nfeval, t_pop *tpa_array, int i_NP)
{
	
  std::cout.precision(16);
  double *p0 = new double[5], f_result=0.0;
  for(int i = 0; i< i_D; i++ )
  {
    p0[i] = t_tmp.fa_vector[i];
  }  
  f_result = c_main(p0,i_D);
  t_tmp.fa_cost[0] = f_result;  
  return(t_tmp);  
}


//------------------Selection criterion of DE--------------------- 
/*
* Description:     Selection criterion of DE. Decide when the trial
*                  vector wins over the target vector
*
* Parameters:      t_trial        trial vector
*                  t_target       target vector
*
* Return value:    TRUE if trial vector wins, FALSE otherwise
*
*/                                                                  
int left_vector_wins(t_pop t_trial, t_pop t_target)
{
	//---trial wins against target even when cost is equal.-----
	if (t_trial.fa_cost[0] <= t_target.fa_cost[0]) return(TRUE);
	else return(FALSE);
}


