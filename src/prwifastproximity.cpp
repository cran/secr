#include <Rcpp.h>
#include <RcppParallel.h>
#include "secr.h"

using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================
// 2019-09-05, 
// 2020-04-05 fixed bug in indiv option
// detector type count
// one occasion, binomial size from integer Tsk
// no deaths, no groups

struct fasthistories : public Worker {
    
    // input data
    const int   mm;
    const int   nc;
    const int   cc; // number of parameter combinations
    const int   grain;
    const int   binomN;
    const bool  indiv;
    const RMatrix<int>    w;          // n x k
    const RMatrix<int>    ki;         // n x k
    const RVector<double> gk; 
    const RVector<double> hk; 
    const RVector<double> density;    // n
    const RVector<int>    PIA;
    const RVector<int>    Tsk;        // k
    const RMatrix<int>    mbool;      // appears cannot use RMatrix<bool>
    
    // internal work arrays
    std::vector<double> pm0base;
    std::vector<double> pm0kbase;
    
    int kk;
    
    // output likelihoods
    RVector<double> output;
    
    // Constructor to initialize an instance of Somehistories 
    // The RMatrix class can be automatically converted to from the Rcpp matrix type
    fasthistories(
        const int mm, 
        const int nc, 
        const int cc,
        const int grain,                    
        const int binomN,
        const bool indiv,
        const IntegerMatrix w,
        const IntegerMatrix ki,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericVector density,
        const IntegerVector PIA,
        const IntegerVector Tsk,
        const LogicalMatrix mbool,
        NumericVector output)    
        : 
        mm(mm), nc(nc), cc(cc), grain(grain), 
        binomN(binomN), indiv(indiv),
        w(w), ki(ki), gk(gk), hk(hk), density(density), PIA(PIA), Tsk(Tsk), mbool(mbool),
        output(output) {
        kk = Tsk.size();   // assuming single occasion
        
        // initialise base value of pm arrays (for n = 0)
        pm0base.reserve(mm);   // for std::vector
        pm0kbase.reserve(kk*mm);
        pr0(0, pm0base, pm0kbase);
    }
    //==============================================================================
    
    void pr0 (const int n, std::vector<double> &pm0n, std::vector<double> &pm0kn) { 
            int c, k, m, w3;
        for (m=0; m<mm; m++) pm0n[m] = 1.0;
        for (k=0; k<kk; k++) {
            w3 =  i3(n, 0, k, nc, 1);    // allow individual or trap covariate 2019-12-06
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (m=0; m<mm; m++) {
                    if (binomN==0)
                        pm0kn[kk * m + k] = gpois (0, Tsk[k]* hk[i3(c, k, m, cc, kk)], 0);
                    else
                        pm0kn[kk * m + k] = gbinom (0, Tsk[k], gk[i3(c, k, m, cc, kk)], 0);
                    pm0n[m] *= pm0kn[kk * m + k];
                }
            }
        }
    }
    //==============================================================================
    
    void prwL (const int n, std::vector<double> &pm) {
        int c, i, k, m, w3;
        
        std::vector<double> pm0(mm);
        std::vector<double> pm0k(kk*mm);
        double pm0ktmp;
        bool base = (n==0) || !indiv;
        
        if (!base) { // (n>0) && indiv) {
            pr0(n, pm0, pm0k); // update for this animal - expt 2019-12-06, 2020-04-04
        }

        if (base) 
            for (m=0; m<mm; m++) pm[m] = pm0base[m];
        else 
            for (m=0; m<mm; m++) pm[m] = pm0[m];  
        for (i=0; i<kk; i++) {
            k = ki(n,i);
            if (k<0) break;    // no more sites
            w3 =  i3(n, 0, k, nc, 1);
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (m=0; m<mm; m++) {
                    if (mbool(n,m)) {
                        if (base) 
                            pm0ktmp = pm0kbase[kk * m + k]; 
                        else 
                            pm0ktmp = pm0k[kk * m + k];
                        if (binomN==0)
                            pm[m] *= gpois (w(n,i), Tsk[k]*hk[i3(c, k, m, cc, kk)], 0) / pm0ktmp;
                        else
                            pm[m] *= gbinom (w(n,i), Tsk[k], gk[i3(c, k, m, cc, kk)], 0) / pm0ktmp;
                    }
                    else {
                        pm[m] = 0.0; 
                    }
                }
            }
        }
    }
    //==============================================================================
    
    double onehistory (int n) {
        std::vector<double> pm(mm);
        prwL (n, pm);           
        for (int m=0; m<mm; m++) {
            pm[m] *= density[m];
        }
        double value = std::accumulate(pm.begin(), pm.end(), 0.0); // may be zero
        return value;
    }
    //==============================================================================
    
    // function call operator that works for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {        
        for (std::size_t n = begin; n < end; n++) {
            output[n] = onehistory (n);
        }
    }
    //==============================================================================
};

// [[Rcpp::export]]
NumericVector fasthistoriescpp (
        const int mm, 
        const int nc, 
        const int cc, 
        const int grain, 
        const int binomN,
        const bool indiv,
        const IntegerMatrix w,
        const IntegerMatrix ki,
        const NumericVector gk, 
        const NumericVector hk, 
        const NumericVector density,
        const IntegerVector PIA, 
        const IntegerVector Tsk,
        const LogicalMatrix mbool) {
    
    NumericVector output(nc); 

    // Construct and initialise
    fasthistories somehist (mm, nc, cc, grain, binomN, indiv,
                            w, ki, gk, hk,
                            density, PIA, Tsk, mbool, 
                            //pm0base, pm0kbase, 
                            output); 
    
    if (grain>0) {
        // Run operator() on multiple threads
        parallelFor(0, nc, somehist, grain);
    }
    else {
        // for debugging avoid multithreading and allow R calls e.g. Rprintf
        somehist.operator()(0,nc);    
    }
    
    // Return consolidated result
    return output;
}
//==============================================================================
