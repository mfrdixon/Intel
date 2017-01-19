/*
* Declarations for error function
* and option class for storing the option data
*
* Author: Matthew Dixon and Xin Zhang
*
* Date: Nov.30 2016
*/

#ifndef __HESTONFOURIERCOSINE_HPP__
#define __HESTONFOURIERCOSINE_HPP__

double Error_Function(
        double* p0,     		//parameters
        double s,       		//option underlying
        double* K,      		//option strike
        double* T,      		//option maturity
        char* Types,    		//option type {'C','P'}
        double* OP,     		//observed prices
        double* W,      		//weights
        int n,          		//chain size
        int nInt,       		//number of terms in fourier-cosine series
        bool bIncorporateNLContraint);  // encode non-linear constraint

class Option{
    double strike;
    double maturity;
    double bid;
    double ask;
    double underlying;
    char type;
public:
    Option(double strike, double maturity, double bid, double ask, double underlying, char type): strike(strike), maturity(maturity), bid(bid), ask(ask), underlying(underlying), type(type){;};
    double getMidPrice() { return (bid+ask)/2.0;};
    double getWeight() { return (((ask-bid) == 0.0) ? 0.001 : (1.0/(ask-bid))); };
    double getStrike(){ return strike;};
    double getMaturity(){ return maturity;};
    double getUnderlying(){ return underlying;};
    double getType(){ return type;};
};
#endif
