#include <Rcpp.h>
using namespace Rcpp;

#include "secr.h"

// [[Rcpp::export]]
List getdenomcpp (int fn, 
                  NumericVector miscparm, 
                  NumericMatrix mask, 
                  int mm,
                  double sigma, 
                  double z) 

{
  double sigma2, lam, dsq;
  double d = 1.0;
  double expzj = 1.0;
  double denomtot = 0.0;
  double scale = 1.0;
  std::vector<double> invdenom(mm);
  
  int m, j;
  int checkinterval = 100;
  //---------------------------------------------------------
  // normalization 
  if (((fn == 4) || ((fn >= 14) && (fn <= 19))) && (fabs(miscparm[0]) > 0.5)) {
    sigma2 =  sigma * sigma;
    for (m=0; m<mm; m++) {
      invdenom[m] = 0;
      for (j=0; j<mm; j++) {
        if ((fabs(miscparm[1]) > 0.5))
          expzj = mask[3 * mm + j];  // exp(covariate val) at mask point j 
        dsq = d2cpp(j, m, mask, mask);
        if (fn != 14) d = sqrt(dsq);
        if (fn == 14) lam = exp(-dsq / 2 / sigma2);
        else if (fn == 4) lam = (double) (d <= sigma);
        else if (fn == 15) lam = 1 - exp(- pow(d /sigma , - z));
        else if (fn == 16) lam = exp(-d / sigma);
        else if (fn == 17) lam = exp(-(d-z)*(d-z) / 2 / sigma2);
        else if (fn == 18) lam = R::pgamma(d,z,sigma/z,0,0);
        else if (fn == 19) lam = exp(- pow(d /sigma , z));
        else stop("unrecognised fn");
        invdenom[m] += lam * expzj;
      }
      denomtot += invdenom[m];
      if (m % checkinterval == 0)
        R_CheckUserInterrupt();
    }
    // optional scaling 
    if (fabs(miscparm[2]) > 0.5) {
      scale = denomtot / mm;
      for (m = 0; m < mm; m++) {
        if (invdenom[m] > 0)
          invdenom[m] = scale / invdenom[m];
      }
    }
  }	
  return (List::create(Named("invdenom") = wrap(invdenom), 
                       Named("scale") = scale));
}
//=============================================================
