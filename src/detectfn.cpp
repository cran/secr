#include <Rcpp.h>
#include "secr.h"
using namespace Rcpp;

//     Detection functions
//     2019-07-29 C++ and consolidate
//     2019-10-13 eliminate unused functions
// 
// The function pfn() selects gfnr on the fly.
// pfn() is used in pdotpoint, naiveRPSV, simdetect, trappingsingle, trappingmulti,
// trappingproximity, trappingcount
// 
// All g--- functions return the detection probability g(x). For count detectors this
// must be transformed to the cumulative hazard -log(1-g(x)) (see function 'hazard' 
// in utils.c).
// 
// The gfnr is used 
// 
// (i) for polygon and transect detectors, 
// (ii) for trappingXXX simulations  (via pfn)
// 
// for which the function must be calculated on the fly within the integration algorithm 
// rather than from a lookup precomputed for a fixed set of trap-mask distances.
// 
// For each form there is a corresponding selection function:
//  getzfnr   polygon.cpp, telemetry.cpp
//  getgfnr   pdot.cpp, trapping.cpp (pfn)
//
//                                    getzfnr    getgfnr
//
//                                    fnptr      fnptr  
//                                              
// fn 0  halfnormal                   zhnr      ghnr    
// fn 1  hazard-rate                  zhrr      ghrr    
// fn 2  exponential                  zexr      gexr    
// fn 3  compound halfnormal          zhncr     ghncr   
// fn 4  uniform                      zunr      gunr    
// fn 5  w-exponential                zhfr      ghfr    
// fn 6  annular normal               zhanr     ganr    
// fn 7  cumulative lognormal         zclnr     gclnr   
// fn 8  cumulative gamma             zcgr      gcgr    
// fn 9  binary signal strength       zhsigbinr gsigbinr
// fn 10 signal strength              zhsigr    gsigr   
// fn 11 signal strength + ss         zhsigsphr gsigsph
// fn 12 signal strength + noise      --        gsigSNr
// fn 13 signal strength + ss + noise --        gsigsphSNr 
// fn 14 hazard halfnormal            zhhnr     ghhnr      
// fn 15 hazard hazard rate           zhhrr     ghhrr      
// fn 16 hazard exponential           zhexr     ghexr      
// fn 17 hazard annular normal        zhanr     ghanr      
// fn 18 hazard cumulative gamma      zhcgr     ghcgr      
// fn 19 hazard variable power        etc.

//--------------------------------------------------------------------
// define functions with third parameter z

int par3 (int fn) {
    if ((fn==1) || (fn==3) || (fn == 5)  || (fn == 6)  || (fn == 7) || 
	(fn == 8) || (fn==10) || (fn == 11)  || (fn == 12)  || (fn == 13) || 
        (fn == 15) || (fn==17) || (fn == 18))
	return(1);
    else
	return(0);
}

//--------------------------------------------------------------------
double zhnr (const NumericVector param, const double r) {
    return(-log(1-param[0] * exp(- r * r / 2 / param[1] / param[1])));
}
//--------------------------------------------------------------------
double zhrr (const NumericVector param, const double r) {
    return(-log(1-param[0] * (1 - exp(- pow(r / param[1], -param[2])))));
}
//--------------------------------------------------------------------
double zexr (const NumericVector param, const double r) {
    return (-log(1-param[0] * exp(-r / param[1])));
}
//--------------------------------------------------------------------
double zhncr (const NumericVector param, const double r) {
    double temp;
    temp = param[0] * exp(- r * r  / 2 / param[1] / param[1]);
    if (round(param[2]) > 1) temp = 1 - pow(1 - temp, param[2]);
    return (-log(1-temp));
}
//--------------------------------------------------------------------
double zunr (const NumericVector param, const double r) {
    if (r<param[1]) return (-log(1-param[0]));
    else return (0);
}
//--------------------------------------------------------------------
double zhfr (const NumericVector param, const double r) {
    if (r<param[2]) return (param[0]);
    else return (-log(1-param[0] * exp(-(r-param[2]) / param[1])));
}
//--------------------------------------------------------------------
double zanr (const NumericVector param, const double r) {
    return (-log(1-param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
				  param[1] / param[1])));
}
//--------------------------------------------------------------------
double zclnr (const NumericVector param, const double r) {
    double g0, sigma, z, CV2, meanlog, sdlog;
    g0 = param[0];
    sigma = param[1];
    z = param[2];
    CV2 = z*z/sigma/sigma;
    meanlog = log(sigma) - log(1 + CV2)/2;
    sdlog = sqrt(log(1 + CV2));
    return (-log(1-g0 * R::plnorm(r,meanlog,sdlog,0,0))); 
}
//--------------------------------------------------------------------
double zcgr (const NumericVector param, const double r) {
    return (-log(1-param[0] * R::pgamma(r,param[2],param[1]/param[2],0,0))); 
}
//--------------------------------------------------------------------

// binary signal strength - (beta0-c)/sdS, beta1/sdS 
double zsigbinr  (const NumericVector param, const double r) {
    double gam, b0, b1;
    b0 = param[0];
    b1 = param[1];
    gam = -(b0 + b1 * r);
    return (-log(1-R::pnorm(gam,0,1,0,0)));    // upper 
}
//--------------------------------------------------------------------

// signal strength - beta0, beta1, sdS 
double zsigr (const NumericVector param, const double r) {
    double mu, gam;
    double beta0, beta1, sdS, cut;
    beta0 = param[0];
    beta1 = param[1];
    sdS = param[2];
    cut = param[3];
    mu = beta0 + beta1 * r;
    gam = (cut - mu) / sdS;
    return (-log(1-R::pnorm(gam,0,1,0,0)));    // upper 
}
//--------------------------------------------------------------------

// signal strength with spherical spreading - beta0, beta1, sdS 
double zsigsphr (const NumericVector param, const double r) {
    double mu, gam;
    double beta0, beta1, sdS, cut;
    beta0 = param[0];
    beta1 = param[1];
    sdS = param[2];
    cut = param[3];
    mu =  beta0 + beta1 * (r-1) - 10 * log(r*r) / M_LN10;
    gam = (cut - mu) / sdS;
    return (-log(1-R::pnorm(gam,0,1,0,0)));    // upper 
}
//--------------------------------------------------------------------

// hazard halfnormal 
double zhhnr (const NumericVector param, const double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard hazard rate 
double zhhrr (const NumericVector param, const double r) {
    return(param[0] * ( 1 - exp(- pow(r / param[1], -param[2]))));
}
//--------------------------------------------------------------------

// hazard exponential 
double zhexr (const NumericVector param, const double r) {
    return (param[0] * exp(-r / param[1]));
}
//--------------------------------------------------------------------

// hazard annular normal 
double zhanr (const NumericVector param, const double r) {
    return (param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
				   param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard cumulative gamma 
double zhcgr (const NumericVector param, const double r) {
    return ((1 - exp( - param[0] * exp(-r / param[1]))));
}

// hazard variable power 
double zhvpr (const NumericVector param, const double r) {
    return (param[0] * exp( - pow(r / param[1], param[2]))) ;
}
//--------------------------------------------------------------------

// hazard halfnormal 
double zhhnrC (const std::vector<double> param, const double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard hazard rate 
double zhhrrC (const std::vector<double> param, const double r) {
    return(param[0] * ( 1 - exp(- pow(r / param[1], -param[2]))));
}
//--------------------------------------------------------------------

// hazard exponential 
double zhexrC (const std::vector<double> param, const double r) {
    return (param[0] * exp(-r / param[1]));
}
//--------------------------------------------------------------------

// hazard annular normal 
double zhanrC (const std::vector<double> param, const double r) {
    return (param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
            param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard cumulative gamma 
double zhcgrC (const std::vector<double> param, const double r) {
    return ((1 - exp( - param[0] * exp(-r / param[1]))));
}

// hazard variable power 
double zhvprC (const std::vector<double> param, const double r) {
    return (param[0] * exp( - pow(r / param[1], param[2]))) ;
}


//====================================================================

double mufn (
    const int k,
    const int m,
    const double b0,
    const double b1,
    const NumericMatrix &A1,
    const NumericMatrix &A2,
    const int spherical)
  
   // Return predicted signal strength at m for source at point k,
   // given strength at source of b0 dB and attenuation of b1 dB/m.
   // Spherical spreading is included if spherical > 0
   // Coordinates of points are in A1 and A2 which have respectively
   // A1rows and A2rows

{
    double d2val;
    d2val = d2cpp(k,m, A1, A2);
    if (spherical <= 0)
	return (b0 + b1 * sqrt(d2val));
    else {
	if (d2val>1) {
	    return (b0 - 10 * log ( d2val ) / 2.302585 + b1 * (sqrt(d2val)-1)); 
	}
	else
	    return (b0);
    }

}
//==============================================================================

double mufnL (
    const int k,
    const int m,
    const double b0,
    const double b1,
    const NumericMatrix &dist2,
    const int spherical)
   // Return predicted signal strength at k for source at point m,
   // given strength at source of b0 dB and attenuation of b1 dB/m.
   // Spherical spreading is included if spherical > 0
   // Uses distance lookup in dist2
   
{
  double d2val;
  d2val = dist2(k,m);
  if (spherical <= 0)
    return (b0 + b1 * sqrt(d2val));
  else {
    if (d2val>1) {
      return (b0 - 10 * log ( d2val ) / 2.302585 + b1 * (sqrt(d2val)-1)); 
    }
    else
      return (b0);
  }
  
}
//==============================================================================

double ghnr (const NumericVector param, const double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
//--------------------------------------------------------------------
double ghrr (const NumericVector param, const double r) {
    return(param[0] * (1 - exp(- pow(r / param[1], -param[2]))));
}
//--------------------------------------------------------------------
double gexr (const NumericVector param, const double r) {
    return (param[0] * exp(-r / param[1]));
}
//--------------------------------------------------------------------
double ghncr (const NumericVector param, const double r) {
    double temp;
    temp = param[0] * exp(- r * r  / 2 / param[1] / param[1]);
    if (round(param[2]) > 1) temp = 1 - pow(1 - temp, param[2]);
    return (temp);
}
//--------------------------------------------------------------------
double gunr (const NumericVector param, const double r) {
    if (r<param[1]) return (param[0]);
    else return (0);
}
//--------------------------------------------------------------------
double ghfr (const NumericVector param, const double r) {
    if (r<param[2]) return (param[0]);
    else return (param[0] * exp(-(r-param[2]) / param[1]));
}
//--------------------------------------------------------------------
double ganr (const NumericVector param, const double r) {
    return (param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
        param[1] / param[1]));
}
//--------------------------------------------------------------------
double gclnr (const NumericVector param, const double r) {
    double g0, sigma, z, CV2, meanlog, sdlog;
    g0 = param[0];
    sigma = param[1];
    z = param[2];
    CV2 = z*z/sigma/sigma;
    meanlog = log(sigma) - log(1 + CV2)/2;
    sdlog = sqrt(log(1 + CV2));
    return g0 * R::plnorm(r,meanlog,sdlog,0,0); 
}
//--------------------------------------------------------------------
double gcgr (const NumericVector param, const double r) {
    return param[0] * R::pgamma(r,param[2],param[1]/param[2],0,0); 
}
//--------------------------------------------------------------------

// binary signal strength - (beta0-c)/sdS, beta1/sdS 
double gsigbinr  (const NumericVector param, const double r) {
    double gam, b0, b1;
    b0 = param[0];
    b1 = param[1];
    gam = -(b0 + b1 * r);
    return (R::pnorm(gam,0,1,0,0));    // upper 
}
//--------------------------------------------------------------------

// signal strength - beta0, beta1, sdS 
double gsigr (const NumericVector param, const double r) {
    double mu, gam;
    double beta0, beta1, sdS, cut;
    beta0 = param[0];
    beta1 = param[1];
    sdS = param[2];
    cut = param[3];
    mu = beta0 + beta1 * r;
    gam = (cut - mu) / sdS;
    return (R::pnorm(gam,0,1,0,0));    // upper 
}
//--------------------------------------------------------------------

// signal strength with spherical spreading - beta0, beta1, sdS 
double gsigsphr (const NumericVector param, const double r) {
    double mu, gam;
    double beta0, beta1, sdS, cut;
    beta0 = param[0];
    beta1 = param[1];
    sdS = param[2];
    cut = param[3];
    mu =  beta0 + beta1 * (r-1) - 10 * log(r*r) / M_LN10;
    gam = (cut - mu) / sdS;
    return (R::pnorm(gam,0,1,0,0));    // upper 
}
//--------------------------------------------------------------------

// hazard halfnormal 
double ghhnr (const NumericVector param, const double r) {
    // Rprintf("%6.3f %6.3f %6.3f \n", r, param[0], param[1]); 
    return(1 - exp( - param[0] * exp(- r * r / 2 / param[1] / param[1])));
}
//--------------------------------------------------------------------

// hazard hazard rate 
double ghhrr (const NumericVector param, const double r) {
    return(1 - exp( - param[0] * ( 1 - exp(- pow(r / param[1], -param[2])))));
}
//--------------------------------------------------------------------

// hazard exponential 
double ghexr (const NumericVector param, const double r) {
    return (1 - exp( - param[0] * exp(-r / param[1])));
}
//--------------------------------------------------------------------

// hazard annular normal 
double ghanr (const NumericVector param, const double r) {
    return (1 - exp( - param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
				   param[1] / param[1])));
}
//--------------------------------------------------------------------

// hazard cumulative gamma 
double ghcgr (const NumericVector param, const double r) {
    return (1 - exp( - (1 - exp( - param[0] * exp(-r / param[1])))));
}
//--------------------------------------------------------------------

// hazard variable power 
double ghvpr (const NumericVector param, const double r) {
    return(1 - exp( - param[0] * exp(- pow(r / param[1], param[2]))));
}
//====================================================================

double pfn (
    const int fn,
    const double d2val,
    const NumericVector &gsb,
    const NumericVector &miscparm,
    const double w2)
{
  double p = -1;
  fnptr gfnr;
  NumericVector tmp(4);
  
  if (d2val > w2) 
    p = 0;
  else {
    gfnr = getgfnr (fn);
    tmp(0) = gsb(0);
    tmp(1) = gsb(1);
    tmp(2) = gsb(2);
    tmp(3) = miscparm[0];
    p = gfnr (tmp, sqrt(d2val));
  }
  return (p);
}

fnptr getgfnr (const int fn) 
{
  if (fn == 0)
    return(ghnr);
  else if (fn == 1)
    return(ghrr);
  else if (fn == 2)
    return(gexr);
  else if (fn == 3)
    return(ghncr);
  else if (fn == 4)
    return(gunr);
  else if (fn == 5)
    return(ghfr);
  else if (fn == 6)
    return(ganr);
  else if (fn == 7)
    return(gclnr);
  else if (fn == 8)
    return(gcgr);
  else if (fn == 9)
    return(gsigbinr);
  else if (fn == 10)
    return(gsigr);
  else if (fn == 11)
    return(gsigsphr);
  else if (fn == 12)
    return(gsigsphr);
  else if (fn == 14)
    return(ghhnr);
  else if (fn == 15)
    return(ghhrr);
  else if (fn == 16)
    return(ghexr);
  else if (fn == 17)
    return(ghanr);
  else if (fn == 18)
    return(ghcgr);
  else if (fn == 19)
    return(ghvpr);
  else // stop("unknown or invalid detection function");
  return(ghnr);
}
//--------------------------------------------------------------------

// For polygon integral 
fnptr getzfnr (const int fn) 
{
  if (fn == 0)
    return(zhnr);
  else if (fn == 1)
    return(zhrr);
  else if (fn == 2)
    return(zexr);
  else if (fn == 3)
    return(zhncr);
  else if (fn == 4)
    return(zunr);
  else if (fn == 5)
    return(zhfr);
  else if (fn == 6)
    return(zanr);
  else if (fn == 7)
    return(zclnr);
  else if (fn == 8)
    return(zcgr);
  else if (fn == 9)
    return(zsigbinr);
  else if (fn == 10)
    return(zsigr);
  else if (fn == 11)
    return(zsigsphr);
  else if (fn == 12)
    return(zsigsphr);
  else if (fn == 14)
    return(zhhnr);
  else if (fn == 15)
    return(zhhrr);
  else if (fn == 16)
    return(zhexr);
  else if (fn == 17)
    return(zhanr);
  else if (fn == 18)
    return(zhcgr);
  else if (fn == 19)
    return(zhvpr);
  else (stop("unknown or invalid detection function"));
  return(ghnr);
}
//--------------------------------------------------------------------

// For polygon integral 
fnptrC getzfnrC (const int fn) 
{
    if (fn == 14)
        return(zhhnrC);
    else if (fn == 15)
        return(zhhrrC);
    else if (fn == 16)
        return(zhexrC);
    else if (fn == 17)
        return(zhanrC);
    else if (fn == 18)
        return(zhcgrC);
    else if (fn == 19)
        return(zhvprC);
    else (stop("unknown or invalid detection function"));
    return(zhhnrC);
}
//--------------------------------------------------------------------
