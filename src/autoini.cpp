#include <Rcpp.h>
#include "secr.h"     

using namespace std;
using namespace Rcpp;

//==============================================================================

// 'naive' functions are used to estimate auto initial values
// these use only the halfnormal detection function 4/5/08

// [[Rcpp::export]]
double naivedcpp (
    const double sigma,          // Parameter : detection scale 
    const IntegerVector &wt,      // integer trap weights 
    const NumericMatrix &traps,   // x,y locations of traps (first x, then y)  
    const NumericMatrix &animals, // x,y locations of animals (first x, then y) 
    const int    fn              // code 0 = halfnormal ONLY
)
{
  int    kk;      // number of traps
  int    nc;       // number of animals
  
  double truncate2 = (2.45 * sigma) * (2.45 * sigma);
  double sump  = 0;
  double sumdp = 0;
  double x,y;
  double dij, d21, d22, p1p2;
  int i,j,n;
  
  if (fn != 0)
    stop ("invalid detection function in external function naivedcpp");
  
  kk = traps.nrow();
  nc = animals.nrow();
  
  for (n=0; n<nc; n++)
  {
    x = animals[n];
    y = animals[n + nc];
    
    for (i=0; i<kk; i++) {
      if (wt[i] > 0) {
        for (j=0; j<(i-1); j++) {
          dij = (traps[i] - traps[j]) * (traps[i] - traps[j]) +
            (traps[i+kk] - traps[j+kk]) * (traps[i+kk] - traps[j+kk]);
          d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+kk] - y) * (traps[i+kk] - y);
          d22 = (traps[j] - x) * (traps[j] - x) + (traps[j+kk] - y) * (traps[j+kk] - y);
          
          if ((d21<=truncate2) && (d22<=truncate2)) {
            p1p2 = exp(-(d21+d22) / 2 / sigma / sigma);
            if (wt[i] > 1)
              p1p2 = 1 - pow(1-p1p2, wt[i]);   
          }
          else
            p1p2 = 0;
          sump  += p1p2;
          sumdp += p1p2 * sqrt(dij);
        }
      }
    }
    for (i=0; i<kk; i++) {  // diagonal 
        d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+kk] - y) * (traps[i+kk] - y);
        if (d21<=truncate2)                                     // d21=d22
	    sump += exp(-2*d21 /2 / sigma / sigma)/2;
    }
    
  }
  return(sumdp/sump);
}

//==============================================================================

// [[Rcpp::export]]
double naivecap2cpp (
    const int    detect,        // scalar code 0 = multicatch, 1 = proximity, 2 = count
    const int    binomN,
    const double g0,            // Parameter : detection magnitude 
    const double sigma,         // Parameter : detection scale 
    int    ss,            // number of occasions
    const IntegerVector &wt,    // integer trap weights 
    const NumericMatrix &traps, // x,y locations of traps (first x, then y)  
    const NumericMatrix &mask,  // x,y locations of mask points (first x, then y) 
    const int    fn                   
)
{
  int    kk;      // number of traps
  int    mm;      // number of mask points
  kk = traps.nrow();
  mm = mask.nrow();
  
  double product;
  double d2val;
  double pk;
  int m,k;
  double nsum = 0;
  double psum = 0;
  int varying = 0;
  
  if (fn != 0)
    stop ("invalid detection function in naivecap2cpp");
  
  if (detect==2 && binomN>1)
      ss = binomN;
  
  if (!varying) {
    for (m=0; m<mm; m++)
    {
      product = 1.0;
      for (k=0; k<kk; k++)
      {
        d2val = d2cpp(k, m, traps, mask);
        pk = g0 * expmin(-d2val / 2 / sigma / sigma);
        product *= (1 - pk); 
        if (detect >= 1) nsum += pk;
      }
      if (detect == 0) nsum += (1 - product);
      psum += 1 - pow(product, ss);
    }
    if (psum<=0)
      return(0);    
    else
      return(ss * nsum / psum);
  }
  else {
    
    // abandon multicatch for now as requires more complex allowance for trap competition 
    // for multicatch need to accumulate Pr(caught) within each occasion, and 
    //   need entire usage for this 
    
    for (m=0; m<mm; m++) {
      pk = 0.0;
      product = 1.0;
      for (k=0; k<kk; k++)
      {
        if (wt[k] > 0) {
          d2val = d2cpp(k, m, traps, mask);
          pk = g0 * expmin(-d2val / 2 / sigma / sigma);
          nsum += wt[k] * pk;  // wt[k] opportunities, each probability pk
          if (wt[k] > 1)
            product *= pow(1 - pk, wt[k]); 
          else
            product *= (1 - pk); 
        }
      }
      psum += 1 - product;   // Pr capt animal at m over all traps & times 
    }
    if (psum<=0)
      return(0);    // failed 
    else
      return(nsum / psum);
  }
}

/*==============================================================================*/
