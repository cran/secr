#include <Rcpp.h>
#include <RcppParallel.h>
#include "secr.h"
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

//===============================================================================

struct Hckmpoly : public Worker {
  
  // input data
  const int detectfn;
  const int dim;
  const RMatrix<double> gsbval;
  const RVector<int>    cumk;
  const RMatrix<double> traps;
  const RMatrix<double> mask;
  
  // output vector to write to
  RVector<double> H;
  RVector<double> gk;
  RVector<double> hk;
  
  int cc, kk, nk, mm, npar;
  double *ex; 
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  Hckmpoly(const int detectfn,
           const int dim,
           const NumericMatrix gsbval, 
           const IntegerVector cumk, 
           const NumericMatrix traps, 
           const NumericMatrix mask, 
           NumericVector H,
           NumericVector gk,
           NumericVector hk)
    : detectfn(detectfn), dim(dim), gsbval(gsbval), cumk(cumk), 
      traps(traps), mask(mask), H(H), gk(gk), hk(hk) {
    
    cc = gsbval.nrow();
    kk = traps.nrow();
    mm = mask.nrow();
    nk = cumk.size()-1;
    npar = gsbval.ncol();
    ex = (double *) R_alloc(10 + 2 * traps.nrow(), sizeof(double));
    
  }
  
  double integral1DRcpp
    (const int fn, 
     const int m, 
     const int c, 
     const RMatrix<double> &gsbval, 
     const RMatrix<double> &traps,
     const RMatrix<double> &mask, 
     const int n1, 
     const int n2, 
     double ex[])
  {
    double ax=0;
    double bx=0;
    double epsabs = 0.0001;
    double epsrel = 0.0001;
    double result = 0;
    double abserr = 0;
    int neval = 0;
    int ier = 0;
    int limit = 100;
    int lenw = 400;
    int last = 0;
    int iwork[100];
    double work[400];
    int k;
    int ns;
    ns = n2-n1+1;
    
    // 2011-06-21 uniform treated separately 
    if (fn == 4) {
      for (k=n1+1; k<=n2; k++) {  // upper bound is length of this transect 
        bx += SegCircle2(traps(k-1,0), traps(k-1,1), traps(k,0), traps(k,1),
                         mask(m,0), mask(m,1), gsbval(c,1));
      }
      return (bx);
    }
    
    for (k=n1+1; k<=n2; k++) {  // upper bound is length of this transect 
      bx += sqrt( (traps(k,0) - traps(k-1,0)) * (traps(k,0) - traps(k-1,0)) +
        (traps(k,1) - traps(k-1,1)) * (traps(k,1) - traps(k-1,1)) );
    }
    // pass parameters etc. through pointer 
    ex[0] = gsbval(c,0);
    ex[1] = gsbval(c,1);
    ex[2] = gsbval(c,2);
    ex[3] = fn;
    ex[4] = mask(m,0);
    ex[5] = mask(m,1);
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;
    for (k=0; k<ns; k++) {               // pass transect vertices 
      ex[k+10] = traps(k+n1,0);        // x 
      ex[k+ns+10] = traps(k+n1,1);     // y 
    }
    Rdqags(fx1, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
           &limit, &lenw, &last, iwork, work);
    if (ier != 0) Rprintf("ier error code in integral1Dcpp %5d\n", ier);
    return (result);
  }
  
  double integral2DRcpp  (const int fn, 
                          const int m, 
                          const int c, 
                          const RMatrix<double> &gsbval, 
                          const RMatrix<double> &traps,
                          const RMatrix<double> &mask, 
                          const int n1, 
                          const int n2, 
                          double ex[]) {
    double ax=1e20;
    double bx=-1e20;
    double epsabs = 0.0001;
    double epsrel = 0.0001;
    double result = 0;
    double abserr = 0;
    int neval = 0;
    int ier = 0;
    int limit = 100;
    int lenw = 400;
    int last = 0;
    int iwork[100];
    double work[400];
    int k;
    int ns;
    int reportier = 0;
    
    // limits from bounding box of this polygon 
    ns = n2-n1+1;
    for (k=0; k<ns; k++) {
      ax = R::fmin2(ax, traps[k+n1]);
      bx = R::fmax2(bx, traps[k+n1]);
    }
    
    // pass parameters etc. through pointer 
    ex[0] = gsbval(c,0);
    ex[1] = gsbval(c,1);
    ex[2] = gsbval(c,2);
    ex[3] = fn;
    ex[4] = mask(m,0);
    ex[5] = mask(m,1);
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;
    
    // also pass polygon vertices 
    for (k=0; k<ns; k++) {
      ex[k+10] = traps(k+n1,0);        // x 
      ex[k+ns+10] = traps(k+n1,1);       // y 
    }
    Rdqags(fx, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
           &limit, &lenw, &last, iwork, work);
    if ((ier != 0) & (reportier))
      Rprintf("ier error code in integral2D %5d\n", ier);
    return (result);
  }
  
  // par[0] = gsbval[c];
  // par[1] = gsbval[cc + c];
  // par[2] = gsbval[2* cc + c];
  // H = hintegral1(fn, par);         /* 2017-03-22 unbounded integrated hazard */
  // detspec[2+c] = H;               /* passed to prwitransect */
  // for (k=0; k<nk; k++) {               /* over transects */
  //     for (m=0; m<mm; m++) {
  //         gi = i3(c,k,m,cc,nk);
  //         /* 2017-03-22 strictly use hazard form */
  //         hk[gi] = par[0] * integral1D (fn, m, c, gsbval, cc,
  //                                       traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], mm, ex) / H;
  //         gk[gi] = 1 - exp(-hk[gi]);
  //     }
  // }
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (int c=0; c<cc; c++) {
      NumericVector gsb(3);
      for (int i=0; i<npar; i++) gsb[i] = gsbval(c, i);
      if (dim==1)
        H[c] = hintegral1cpp(detectfn, gsb);         // unbounded integrated hazard from radial function
      else
        H[c] = hintegralcpp(detectfn, gsb);         // unbounded integrated hazard from radial function
      for (int k=0; k<nk; k++) {                  // over parts 
        for (std::size_t m = begin; m < end; m++) {
          int gi = i3(c, k, m, cc, nk);
          // strictly use hazard form : expected detections of animals at m 
          // gsb[0] only makes sense here if fn is HHN, HHR, HEX, HAN HCG 
          if (dim==1)
            hk[gi] = gsb[0] * integral1DRcpp (detectfn, m, 0, gsbval, traps, mask, cumk[k], cumk[k+1]-1, ex) / H[c];
          else 
            hk[gi] = gsb[0] * integral2DRcpp (detectfn, m, 0, gsbval, traps, mask, cumk[k], cumk[k+1]-1, ex) / H[c];
          gk[gi] = 1 - exp(-hk[gi]);
        }
      }
    }
  }
};

// [[Rcpp::export]]
List makegkpolygoncpp (const int detectfn, 
                       const int dim,
                       const int grain,
                       const NumericMatrix& gsbval, 
                       const IntegerVector& cumk,
                       const NumericMatrix& traps,
                       const NumericMatrix& mask
) 
{
  NumericVector H(gsbval.nrow()); 
  NumericVector gk(gsbval.nrow() * (cumk.size()-1) * mask.nrow()); 
  NumericVector hk(gsbval.nrow() * (cumk.size()-1) * mask.nrow()); 
  
  Hckmpoly hckm (detectfn, dim, gsbval, cumk, traps, mask, H, gk, hk);
  
  // Cannot use R call Rdqags if multithreading
  // Restrict to single-thread until understand and apply RcppParallel/RcppNumerical for
  // 2-D integration
  // if (grain>0) {
  //    // call it with parallelFor
  //  parallelFor(0, mask.nrow(), hckm, grain);
  // }
  // else {
  hckm.operator()(0,mask.nrow());    // for debugging avoid multithreading to allow R calls
  // }
  return List::create(Named("H") = H, Named("gk") = gk, Named("hk") = hk);
}
//==============================================================================
