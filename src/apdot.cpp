#include <Rcpp.h>
#include "poly.h"
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

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

// // [[Rcpp::export]]
// NumericVector hdotpolycpp (
//         const NumericMatrix xy, 
//         const NumericMatrix traps, 
//         const IntegerVector detect, 
//         const NumericMatrix Tsk, 
//         const IntegerVector markocc, 
//         const int nk, 
//         const int ss, 
//         const IntegerVector kk,
//         const int fn, 
//         const NumericVector gsb,
//         const bool convex)
// {
//     int i,k,s;
//     double hk = 0;
//     double sumhk;
//     double H = 1.0;
//     double Tski = 1.0;
//     bool allsighting = true;
//     NumericMatrix gsbval(1,gsb.size());
//     int nxy = xy.nrow();
//     std::vector<double> value(nxy);
//     std::vector<int> cumk(nk+1);
//     
//     // ex is used only by Rdqags in integral2Dcpp
//     // double *ex;
//     // ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
//     
//     for (s=0; s<ss; s++) {
//         if (markocc[s]>0) allsighting = false;                    // capture occasions 
//     }
//     
//     cumk[0] = 0;
//     for (i=0; i<nk; i++)
//         cumk[i+1] = cumk[i] + kk[i];
//     
//     // 2017-03-22 integrated hazard 
//     for (i=0;i<gsb.size(); i++) gsbval(0,i) = gsb(i);
//     if (anypolygon(detect))          
//         H = hintegral2Ncpp(fn, as<std::vector<double>>(gsb));     
//     else if (anytransect(detect))    
//         H = hintegral1Ncpp(fn, as<std::vector<double>>(gsb));
//     else 
//         stop ("unrecognised detector type in hdotpolycpp");
//     
//     // conversions for RMatrix input to integralxDNRcpp
//     const RcppParallel::RMatrix<double> gsbvalR(gsbval);
//     const RcppParallel::RMatrix<double> trapsR(traps);
//     const RcppParallel::RMatrix<double> xyR(xy);
//     
//     for (i=0; i<nxy; i++) {
//         sumhk = 0.0;
//         for (s=0; s<ss; s++) {
//             if (((markocc[s]>0) || allsighting) && (detect[s]!=13)) {                     
//                 for (k=0; k<nk; k++) {
//                     Tski = Tsk[s * nk + k];
//                     // subset traps rows cumk[k] : (cumk[k+1]-1)
//                     int n1 = cumk[k];
//                     int n2 = cumk[k+1]-1;
//                     
//                     if (Tski > 1e-10) {  
//                         if (anypolygon(detect))
//                             hk = gsb(0) * integral2DNRcpp (fn, i, 0, gsbvalR, trapsR, xyR, n1, n2, convex) / H;
//                         //hk = gsb(0) * integral2Dcpp (fn, i, 0, gsbvalR, trapsR, xyR, n1, n2, ex) / H;
//                         else if (anytransect(detect))
//                             //hk = gsb(0) * integral1DNcpp (fn, i, 0, gsbval, traps, xy, n1, n2) / H;
//                             hk = gsb(0) * integral1DNRcpp (fn, i, 0, gsbvalR, trapsR, xyR, n1, n2) / H;
//                         sumhk += hk * Tski;
//                     }
//                 }
//             }
//         }
//         value[i] = sumhk;
//     }
//     return(wrap(value));
// }

//=============================================================
struct hdotpoly : public Worker {
    
    // input data
    const int detectfn;
    const bool convex;
    const int dim;
    const RVector<double> gsbR;
    const RMatrix<double> gsbvalR;
    const RVector<int> cumk;
    const RVector<int> markocc;
    const RMatrix<double> trapsR;
    const RMatrix<double> xyR;
    const RMatrix<double> TskR;
    RVector<double> hdot;
    
    double H;
    
    // output vector to write to
    
    int kk, nk, mm, npar, ss;
    bool allsighting = true;

    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    hdotpoly(const int detectfn,
             const bool convex,
             const int dim,
             const NumericVector gsb, 
             const NumericMatrix gsbval, 
             const IntegerVector cumk,
             const IntegerVector markocc,
             const NumericMatrix traps, 
             const NumericMatrix xy, 
             const NumericMatrix Tsk, 
             NumericVector hdot)
        : detectfn(detectfn), convex(convex), dim(dim), gsbR(gsb), gsbvalR(gsbval),
          cumk(cumk), markocc(markocc), trapsR(traps), xyR(xy), TskR(Tsk), 
          hdot(hdot) {
        
        nk = cumk.size()-1;
        npar = gsb.size();
        ss = Tsk.ncol();
        
        for (int s=0; s<ss; s++) {
            if (markocc[s]>0) allsighting = false;                    // capture occasions 
        }
        
        if (dim==1)
            H = hintegral1Ncpp(detectfn, as<std::vector<double>>(gsb));
        else
            H = hintegral2Ncpp(detectfn, as<std::vector<double>>(gsb));     
        
    }
    
    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end)
    {
        double hk = 0.0;
        double sumhk;
        double Tski;
        
        for (std::size_t i = begin; i<end; i++) {
            sumhk = 0.0;
            for (int s=0; s<ss; s++) {
                if ((markocc[s]>0) || allsighting) {                     
                    for (int k=0; k<nk; k++) {
                    Tski = TskR(k,s);
                    // subset traps rows cumk[k] : (cumk[k+1]-1)
                    int n1 = cumk[k];
                    int n2 = cumk[k+1]-1;
                    
                    if (Tski > 1e-10) {  
                        if (dim == 1)
                            hk = gsbR[0] * integral1DNRcpp (detectfn, i, 0, gsbvalR, trapsR, xyR, n1, n2) / H;
                        else // dim == 2
                            hk = gsbR[0] * integral2DNRcpp (detectfn, i, 0, gsbvalR, trapsR, xyR, n1, n2, convex) / H;
                        sumhk += hk * Tski;
                    }
                }
                }
            }
            hdot[i] = sumhk;
        }
    }
};

// [[Rcpp::export]]
NumericVector hdotpolycpp2 (
        const NumericMatrix &xy, 
        const NumericMatrix &traps, 
        const NumericMatrix &Tsk, 
        const IntegerVector &markocc, 
        const IntegerVector &cumk, 
        const int &detectfn, 
        const NumericVector &gsb,
        const bool &convex,
        const int &dim,
        const int &grain)
{
    int nxy = xy.nrow();
    NumericMatrix gsbval(1,gsb.size());
    NumericVector hdot(nxy);
    for (int i=0; i<gsb.size(); i++) gsbval(0,i) = gsb(i);

    hdotpoly hpoly (detectfn, convex, dim, gsb, gsbval, cumk, markocc, traps, xy, Tsk, hdot);
    
    if (grain>0) {
        // call it with parallelFor
        parallelFor(0, nxy, hpoly, grain);
    }
    else {
        hpoly.operator()(0,nxy);    // for debugging avoid multithreading to allow R calls
    }
    
    return(wrap(hdot));
}

//=============================================================

