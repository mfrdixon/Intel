#include "svd.h"
//Fitting routine for an expansion with
//nl Legendre polynomials pl, evaluated using the recurrence relation
void fleg(double x, double pl[], int nl)
{
   int j;
   double twox,f2,f1,d;
   pl[1]=1.0;
   pl[2]=1.0-x;
   if(nl>2){
     for (j=3;j<=nl;j++) {
        f1=j-2;
        f2= 2.0*(j-2)+1-x;
        pl[j]=(f2*pl[j-1]-f1*pl[j-2])/(j-1);
     }
   }
}
