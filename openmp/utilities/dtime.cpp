/*
* Definition of the function for measuring 
* benchmark time
*
* author: Matthew Dixon and Xin Zhang
* Date: Nov. 30 2016
*/

#include <sys/time.h>

double dtime()
{
     double tseconds = 0.0;
     struct timeval mytime;
     gettimeofday(&mytime,(struct timezone*)0);
     tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
     return( tseconds );
}
