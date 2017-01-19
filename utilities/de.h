/*
* Definition of macros and structure
* for population member
* 
* Author: Matthew Dixon and Xin Zhang
* Date: Nov.30 2016
*/


//------------------Prevent multiple includes of de.h-------------------------
#ifndef _DE_H
#define _DE_H

//------------------------General constants-----------------------------------
#define MAXDIM   6         // maximum number of dimensions i.e. parameters. 
#define MAXPOP  2000       // number of random vectors to be stored. Watch  
			      // out! gi_D must be <= 33 because w_index writes 
			      // up to location 3*gi_D-1.                       
#define MAXCOST   1        // maximum number of objectives to be minimized
#define MAXCONST 20        // maximum number of constraints

#define VTR     1.0e-10    // value to reach

#define TRUE  1
#define FALSE 0


//----------------------------Typedefs----------------------------------------
typedef struct
//-----------------Definition of population member----------------------------
{
   double fa_vector[MAXDIM];         //parameter vector
   double fa_cost[MAXCOST];          //vector of objectives (costs)
   double fa_constraint[MAXCONST];   //vector of constraints
} t_pop;

#endif // _DE_H