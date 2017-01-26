#include "HestonFourierCosine.h"
#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <omp.h>

#define pi 3.1415926535897932384626433832795

float counter = 0.0;
float totaltime = 0.0;
std::vector<Option> chain;
double* K ,*T,*W,*OP,s;
char* Types;
int num_blocks;

double dtime();

inline std::complex<double> HestonCF(
        double u,
        double T,
        double r,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho
        ) {
        std::complex<double> j1(0, 1);
    	double a                    = lmbda*meanV;
    	double b                    = lmbda;
    	double sigma2               = sigma*sigma;
    	double e                    = pow(rho*sigma*u-b,2);
    	double d                    = sqrt(e+(u*u+u)*sigma2);
    	std::complex<double> g      = (b-j1*rho*sigma*u - d)/(b-j1*rho*sigma*u+d);
    	std::complex<double> ret    = exp(j1*u*r*T);
    	ret                         *= exp((a/sigma2)*((b - rho*j1*sigma*u - d)*T - 2.0*log((1.0-g*exp(-d*T))/(1.0-g))));
    	return ret*exp((v0/sigma2)*(b - rho*j1*sigma*u - d)*(1.0-exp(-d*T))/(1.0-g*exp(-d*T)));
}
inline double xi(
        double k,
        double a,
		double b,
		double c,
        double d) {

    	double ret = 1.0/(1.0+pow(k*pi/(b-a),2))*(cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a)*sin(k*pi*(d-a)/(b-a))*exp(d)-k*pi/(b-a)*sin(k*pi*(c-a)/(b-a))*exp(c));
    	return ret;
}
inline double psi(
        double k,
        double a,
        double b,
        double c,
        double d) {
    double ret=0.0;
    if (k==0)
       ret = d-c;
    else
       ret = (sin(k*pi*(d-a)/(b-a))-sin(k*pi*(c-a)/(b-a)))*(b-a)/(k*pi);
    return ret;
}
inline double HestonCOS(
        double S,
        double K,
        double T,
        double r,
        double sigma,
        double lmbda,
        double meanV,
        double v0,
        double rho,
        char otype,
        int N) {

    std::complex<double> j1(0, 1);
    double sigma2 = sigma*sigma;
    double lmbda2 = lmbda * lmbda;
    double c1 = r*T+(1-exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T;
    double c2 = 1.0/(8.0*lmbda2*lmbda)*(sigma*T*lmbda*exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma
*(1-exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma2+4.0*lmbda2)+sigma2*((meanV-2.0*v0)
*exp(-2.0*lmbda*T)+meanV*(6.0*exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda2*(v0-meanV)*(1-exp(-lmbda*T)));
    double a = c1-12.0*sqrt(fabs(c2));
    double b = c1+12.0*sqrt(fabs(c2));
    double x = log(S/K);
    double U, unit = 0.5;
    // Note that this for loop is independent of the strike
    double ret = 0.0;
    #pragma omp simd reduction(+:ret)
    for (int k=0; k < N; k++)
    {
      	U = 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b));
      	std::complex<double> HCF = HestonCF(k*pi/(b-a),T,r,sigma,lmbda,meanV,v0,rho);
        ret += unit*(HCF*exp(j1*double(k)*pi*(x-a)/(b-a))).real()*U;
      	unit = 1.0;
    }
    return K*exp(-r*T)*ret;
}

double Error_Function(
        const double* p0,     //parameters
        double s,       //option underlying
        double* K,      //option strike
        double* T,      //option maturity
        char* Types,    //option type {'C','P'}
        double* OP,     //observed prices
        double* W,      //weights
        int n,          //chain size
        int nInt,       //number of terms in fourier-cosine series
	bool bIncorporateNLContraint){  // encode non-linear constraint

      	double kappa, theta, sigma, rho, v0, x;
	double* xtemp = new double[n];
     	double r0 = 0.01;
    	kappa = p0[0]; theta =p0[1]; sigma = p0[2]; rho = p0[3]; v0=p0[4];
    	if (bIncorporateNLContraint)
         kappa = (kappa + sigma*sigma)/(2.0*theta);
        double tmp, r = 0 ;

    	#pragma omp parallel for private(tmp) reduction(+ : r)
    	for (int i=0; i<n; i++){
          tmp = W[i]*(HestonCOS(s,K[i],T[i],r0,sigma,kappa,theta,v0,rho,Types[i],nInt) - OP[i]);
          r += tmp*tmp;
     	}
    	double rmse = sqrt(r/n);
    	return rmse;
}


inline int readIntraDataFromFile(const char* fileName, std::vector<Option> & chain)
{
  std::fstream file(fileName, std::ios::in);
  
  if(!file.is_open())
  {
    std::cout << "File not found!\n";
    return 1;
  }

  std::string csvLine;
  int rowIdx = 0;
  
  while( std::getline(file, csvLine) )
  {
		std::istringstream csvStream(csvLine);
		std::string csvCol;
		std::vector<std::string> csvRow;
		
		while( std::getline(csvStream, csvCol, ',') )
		{
			csvRow.push_back(csvCol);
		}	  
		if (rowIdx++ >0)
		{
			int t = std::atoi(csvRow[1].c_str());
			int T = std::atoi(csvRow[3].c_str());
			double tau = difftime (time_t(T), time_t(t))/ (260.0*24.0*3600.0);
			double strike = std::atof(csvRow[4].c_str());
			double bid = std::atof(csvRow[6].c_str());
			double ask = std::atof(csvRow[7].c_str());
			double underlying = std::atof(csvRow[5].c_str());
			char type = csvRow[2].c_str()[0];
			Option o = Option(strike,tau,bid,ask,underlying,type);
			chain.push_back(o);
		}
  }
  return 0;
}

int taking_input()
{
  	double w = 0.0;
	int idx = 0;
	readIntraDataFromFile("../data/ZNGA-Option_sorted.csv", chain);
	num_blocks = chain.size();
	  
	K = new double[num_blocks];
      	T = new double[num_blocks];
      	W = new double[num_blocks];
      	OP = new double[num_blocks];
      	Types = new char[num_blocks];
	  
	for (std::vector<Option>::iterator it=chain.begin(); it!=chain.end();it++)
    	{
       		K[idx] = it->getStrike();
       		T[idx] = it->getMaturity();
       		w+= it->getWeight();
       		OP[idx] = it->getMidPrice();
       		Types[idx] = it->getType();
       		s = it->getUnderlying();
       		idx++;
    	}
	
	idx = 0;
	
	for (std::vector<Option>::iterator it=chain.begin(); it!=chain.end();it++)
		W[idx++] = it->getWeight()/w;
	return 1;
}

void free()
{
	std::cout<<"Total Number of calls In Error Function= "<<counter<<std::endl;
	std::cout<<"Time Per Error Function Call= "<<totaltime/counter<<" sec"<<std::endl;
	std::cout<<"Total Time Elapsed In Error Function = "<<totaltime<<" sec"<<std::endl;	

    	delete [] K;
    	delete [] T;
    	delete [] W;
    	delete [] OP;
    	delete [] Types;
}



double c_main(unsigned n, const double *p0 ,double *grad, void *data )
{
    	std::cout.precision(16);
    	bool bIncorporateNLConstraint = true;
    	int nInt = 256;
    	double rmse;

	double t1, t2;
    	t1 = dtime();
	rmse = Error_Function(p0,s,K,T,Types,OP,W,num_blocks,nInt,bIncorporateNLConstraint);
	t2 = dtime();
	
    	++counter;
	
	totaltime += t2-t1;	

	return rmse;
}
