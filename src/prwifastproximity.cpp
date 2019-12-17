#include <Rcpp.h>
#include <RcppParallel.h>
#include "secr.h"

using namespace Rcpp;
using namespace RcppParallel;

//==============================================================================
// 2019-09-05
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
    
    // workarrays
    RVector<double> pm0;
    RMatrix<double> pm0k;

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
        NumericVector pm0,
        NumericMatrix pm0k,
        NumericVector output)    
        : 
        mm(mm), nc(nc), cc(cc), grain(grain), 
        binomN(binomN), indiv(indiv),
        w(w), ki(ki), gk(gk), hk(hk), density(density), PIA(PIA), Tsk(Tsk), mbool(mbool),
        pm0(pm0), pm0k(pm0k), output(output) {
        kk = Tsk.size();
        pr0(0, pm0, pm0k);
    }
    //==============================================================================
    
    void pr0 (const int n, NumericVector &pm0, NumericMatrix &pm0k) { 
        int c, k, m, w3;
        for (m=0; m<mm; m++) pm0[m] = 1.0;
        for (k=0; k<kk; k++) {
            // w3 =  i3(0, 0, 0, nc, 1);
            w3 =  i3(n, 0, k, nc, 1);    // allow individual or trap covariate 2019-12-06
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (m=0; m<mm; m++) {
                    if (binomN==0)
                        pm0k(k,m) = gpois (0, Tsk[k]* hk[i3(c, k, m, cc, kk)], 0);
                    else
                        pm0k(k,m) = gbinom (0, Tsk[k], gk[i3(c, k, m, cc, kk)], 0);
                    pm0[m] *= pm0k(k,m);
                }
            }
        }
    }
    //==============================================================================
    
    // ugly repeat of function definition 2019-12-06
    void pr0R (const int n, RVector<double> &pm0, RMatrix<double> &pm0k) { 
        int c, k, m, w3;
        for (m=0; m<mm; m++) pm0[m] = 1.0;
        for (k=0; k<kk; k++) {
            // w3 =  i3(0, 0, 0, nc, 1);
            w3 =  i3(n, 0, k, nc, 1);    // allow individual or trap covariate 2019-12-06
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (m=0; m<mm; m++) {
                    if (binomN==0)
                        pm0k(k,m) = gpois (0, Tsk[k]* hk[i3(c, k, m, cc, kk)], 0);
                    else
                        pm0k(k,m) = gbinom (0, Tsk[k], gk[i3(c, k, m, cc, kk)], 0);
                    pm0[m] *= pm0k(k,m);
                }
            }
        }
    }
    //==============================================================================
    
    void prwL (const int n, std::vector<double> &pm) {
        int c, i, k, m, w3;
        if (n>0 && indiv) {
            pr0R(n, pm0, pm0k); // update for this animal - expt 2019-12-06
        }
        for (int m=0; m<mm; m++) pm[m] = pm0[m];  // assumed missed at all sites...
        for (i=0; i<kk; i++) {
            k = ki(n,i);
            if (k<0) break;    // no more sites
            w3 =  i3(n, 0, k, nc, 1);
            c = PIA[w3] - 1;
            if (c >= 0) {    // ignore unset traps
                for (m=0; m<mm; m++) {
                  if (mbool(n,m)) {
                    if (binomN==0)
                      pm[m] *= gpois (w(n,i), Tsk[k]*hk[i3(c, k, m, cc, kk)], 0) / 
                        pm0k(k,m);
                    else
                      pm[m] *= gbinom (w(n,i), Tsk[k], gk[i3(c, k, m, cc, kk)], 0) / 
                        pm0k(k,m);
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
    NumericVector pm0(mm); 
    NumericMatrix pm0k(Tsk.size(),mm);

    // Construct and initialise
    fasthistories somehist (mm, nc, cc, grain, binomN, indiv,
                            w, ki, gk, hk,
                            density, PIA, Tsk, mbool, pm0, pm0k, output); 
    
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
