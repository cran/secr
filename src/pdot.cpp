#include <Rcpp.h>
using namespace Rcpp;
#include "secr.h"

// [[Rcpp::export]]
NumericVector pdotpointcpp (
    const NumericMatrix xy, 
    const NumericMatrix traps, 
    const NumericMatrix dist2,                   
    const IntegerVector detect, 
    const NumericMatrix Tsk, 
    const IntegerVector markocc, 
    const int fn, 
    const NumericVector gsb, 
    const NumericVector miscparm,
    const double w2, 
    const IntegerVector binomN)
{
  int i,k,s;
  double tempval;
  double p;
  double Tski = 1.0;
  int ss, kk, nxy;
  ss = Tsk.ncol();
  kk = traps.nrow();
  nxy = xy.nrow();
  std::vector<double> value(nxy);
  
  bool allsighting = true;

  if (anypolygon(detect) || anytransect(detect)) {
    stop("pdotpoint not for polygon or transect detectors");
  }
  
  for (s=0; s<ss; s++) {
    if (markocc[s]>0) allsighting = false;                    /* capture occasions */
  }
  
  if (fn>19) stop("pdotpointcpp requires detectfn < 20");
  for (i=0; i<nxy; i++) {
    tempval = 1;
    for (s=0; s<ss; s++) {
      if (((markocc[s]>0) || allsighting) && (detect[s]!=13)) {                     
        for (k=0; k<kk; k++) {
          Tski = Tsk(k,s);
          if (Tski > 1e-10) {
            p = pfn(fn, dist2(k,i), gsb, miscparm, w2);
            /* counts */
            if (detect[s] == 2) {    
              if (binomN[s] == 0)
                p = 1 - countp(0, 0, Tski * hazard(p));
              else if (binomN[s] == 1)
                p = 1 - countp(0, round(Tski), p);
              else {
                if (fabs(Tski-1) > 1e-10)
                  p = 1 - pow(1-p, Tski);
                p = 1 - countp(0, binomN[s], p);
              }
            }
            else {                
              if (fabs(Tski-1) > 1e-10)
                p = 1 - pow(1-p, Tski);
            }
            tempval *= 1 - p;
          }
        }
      }
    }
    value[i] = 1 - tempval;
  }
  return(wrap(value));
}

//=============================================================

// [[Rcpp::export]]
NumericVector hdotpolycpp (
    const NumericMatrix xy, 
    const NumericMatrix traps, 
    const IntegerVector detect, 
    const NumericMatrix Tsk, 
    const IntegerVector markocc, 
    const int nk, 
    const int ss, 
    const IntegerVector kk,
    const int fn, 
    const NumericVector gsb)
{
  int i,k,s;
  double hk = 0;
  double sumhk;
  double H = 1.0;
  double *ex;
  double Tski = 1.0;
  ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
  bool allsighting = true;
  NumericMatrix gsbval(1,gsb.size());
  int nxy;
  nxy = xy.nrow();
  std::vector<double> value(nxy);
  std::vector<int> cumk(nk+1);
  for (s=0; s<ss; s++) {
    if (markocc[s]>0) allsighting = false;                    // capture occasions 
  }
  
  cumk[0] = 0;
  for (i=0; i<nk; i++)
    cumk[i+1] = cumk[i] + kk[i];
  
  // 2017-03-22 integrated hazard 
  for (i=0;i<gsb.size(); i++) gsbval(0,i) = gsb(i);
  if (anypolygon(detect))          
    H = hintegralcpp(fn, gsb);     
  else if (anytransect(detect))    
    H = hintegral1cpp(fn, gsb);
  else 
    stop ("unrecognised detector type in hdotpolycpp");
  
  
  for (i=0; i<nxy; i++) {
    sumhk = 0.0;
    for (s=0; s<ss; s++) {
      if (((markocc[s]>0) || allsighting) && (detect[s]!=13)) {                     
        for (k=0; k<nk; k++) {
          Tski = Tsk[s * nk + k];
          if (Tski > 1e-10) {  
            if (anypolygon(detect))
              hk = gsb(0) * integral2Dcpp (fn, i, 0, gsbval, traps, xy, cumk[k], cumk[k+1]-1, ex) / H;
            else if (anytransect(detect))
              hk = gsb(0) * integral1Dcpp (fn, i, 0, gsbval, traps, xy, cumk[k], cumk[k+1]-1, ex) / H;
            sumhk += hk * Tski;
          }
        }
      }
    }
    value[i] = sumhk;
  }
  return(wrap(value));
}

//=============================================================

