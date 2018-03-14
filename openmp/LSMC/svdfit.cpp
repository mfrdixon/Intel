#include "svd.h"

void svdfit(double x[], double y[], double sig[], int ndata, double a[], int ma,
	double **u, double **v, double w[], double *chisq,
	void (*funcs)(double, double [], int))
{
	void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
	void svdcmp(double **a, int m, int n, double w[], double **v);
	int j,i;
	double wmax,tmp,thresh,sum,*b,*afunc;

	b=dvector(1,ndata);
	afunc=dvector(1,ma);
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		tmp=1.0/(sig[i] + 0.001);
		for (j=1;j<=ma;j++) {
                   u[i][j]=afunc[j]*tmp;
                }
		b[i]=y[i]*tmp;
	}
	svdcmp(u,ndata,ma,w,v);
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=tol*wmax;
	for (j=1;j<=ma;j++)
		if (w[j] < thresh) w[j]=0.0;
	svbksb(u,w,v,ndata,ma,b,a);
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) {
                    sum += a[j]*afunc[j]; 
                }
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
	free_dvector(afunc,1,ma);
	free_dvector(b,1,ndata);
}
