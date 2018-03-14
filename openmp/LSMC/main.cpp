#include "svd.h"

double call_price(const double& S, const double& K, const double& r, const double& v, const double& T);

int main()
{
    double *rnd, chisq, *x, *y, *c, *afunc, *sig, *a, *w, **cvm, **u, **v, **s, **h;

    int i;
    // Model Parameters
    double const s0 = 52.0;  // initial stock level
    double const K = 40.0;   // strike price
    double const T = 1.0;    // time-to-maturity
    double const r = 0.06;   // short rate
    double const sigma = 0.2;// volatility

    // Simulation Parameters
    int const num_steps = 2;  // number of time steps
    int num_points = 10;     // number of simulations
    int const num_coeffs = 3; /* Number of basis functions  */
    double const dt = T / num_steps;
    double const df = exp(-r * dt);
    double const noise = 0.01; /* Estimate of std deviation of measurement error */


    x   = dvector(1, num_points);
    y   = dvector(1, num_points);
    c   = dvector(1, num_points);
    sig = dvector(1, num_points);
    a   = dvector(1, num_coeffs);
    w   = dvector(1, num_coeffs);
    afunc=dvector(1,num_coeffs);
    cvm = dmatrix(1, num_coeffs, 1, num_coeffs);
    u   = dmatrix(1, num_points, 1, num_coeffs);
    v   = dmatrix(1, num_coeffs, 1, num_coeffs);
    s   = dmatrix(1, num_points, 1, num_steps+1);
    h   = dmatrix(1, num_points, 1, num_steps);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31,
                             0x9908b0df, 11,
                             0xffffffff, 7,
                             0x9d2c5680, 15,
                             0xefc60000, 18, 1812433253> generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);
    // Stock Price Paths
    for (int i=1; i<= num_points;i++){
      s[i][1] = s0;
      for (int j=2; j<=(num_steps+1);j++){
         s[i][j] = s[i][j-1]*(1.0+r*dt + sqrt(dt)*sigma*distribution(generator));
         h[i][j] = DMAX(s[i][j]-K,0.0);
      }
      y[i] = h[i][num_steps+1]; // payoff function
      //std::cout<<y[i]<<std::endl;
      sig[i]= y[i]*noise;
    }
    // American Option Valuation by Backwards Induction
    for (int t=num_steps; t>1;t--) {
       for (int i=1; i<= num_points;i++){
         x[i] = s[i][t];
         y[i] *=df;
         //std::cout<<x[i]<<","<<y[i]<<std::endl;
       }
       // Calculate best-fit polynomial coefficients.
       svdfit(x, y, sig, num_points, a, num_coeffs, u, v, w, &chisq, fleg);
       /*
       * Calculate covariance matrix to get a chi-square
       * estimate of goodness of fit
       */
       // diagnostics
       /*svdvar(v, num_coeffs, w, cvm);
       
       std::cout<<"Coefficients for polynomial of order"<< num_coeffs<<std::endl;
       for (i = 1; i <= num_coeffs; i++) {
          std::cout<<a[i]<<"  +-" << sqrt(cvm[i][i])<<std::endl;
       }
       std::cout<<"Chi-squared "<<  chisq <<std::endl;
       */ 
 
       for (i=1;i<=num_points;i++){
         c[i] =0.0;
         fleg(x[i],afunc,num_coeffs); 
         for (int j=1;j<= num_coeffs;j++){ 
           c[i]+= a[j]*afunc[j]; //lagval
         }
         //std::cout<<c[i]<<std::endl;
         // exercise decision
         if (h[i][t]>c[i])
            y[i] = h[i][t];
       }
    }
    double V0=0.0; 
    for (i=1;i<=num_points;i++)
      V0 += y[i]; 
    V0 *= df/num_points;
    std::cout<<"American call option value "<< V0 <<std::endl;

    std::cout<<"European call option value "<<call_price(s0,K,r,sigma,T)<<std::endl;


    /*https://math.stackexchange.com/questions/1939359/least-square-approximation-using-laguerre-polynomials
    for (int i=1; i<=num_points; i++)
    {
        x[i] = i-1;
        y[i] = 9.00/8.0 + (3.0/16.0)*x[i]*x[i];
        sig[i] = noise * y[i];
    }*/
    
    free_dmatrix(s,   1, num_points, 1, num_steps);
    free_dmatrix(h,   1, num_points, 1, num_steps);
    free_dmatrix(v,   1, num_coeffs, 1, num_coeffs);
    free_dmatrix(u,   1, num_points, 1, num_coeffs);
    free_dmatrix(cvm, 1, num_coeffs, 1, num_coeffs);
    free_dvector(w,   1, num_coeffs);
    free_dvector(a,   1, num_coeffs);
    free_dvector(sig, 1, num_points);
    free_dvector(c,   1, num_points);
    free_dvector(y,   1, num_points);
    free_dvector(x,   1, num_points);
    free_dvector(afunc,1,num_coeffs);

    return EXIT_SUCCESS;
}
