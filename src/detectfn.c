/*
    Detection functions
    2011-09-30
    2012-02-01 modified for single mufn
    2014-08-27 mufn moved here from utils.c

Three forms are provided:
hfn  -- distance is an argument (r)
gfn  -- distance is computed from trap amd mask indices and coordinate objects
gfnL -- distance is computed from trap amd mask indices and lookup table

*/
/*--------------------------------------------------------------------*/

/*
fn 0  halfnormal
fn 1  hazard-rate
fn 2  exponential
fn 3  compound halfnormal
fn 4  uniform
fn 5  w-exponential
fn 6  annular normal
fn 7  cumulative lognormal
fn 8  cumulative gamma
fn 9  binary signal strength
fn 10 signal strength
fn 11 signal strength with spherical spreading
fn 12 signal strength + noise
fn 13 signal strength with spherical spreading + noise
fn 14 hazard halfnormal
fn 15 hazard hazard rate
fn 16 hazard exponential
fn 17 hazard annular normal
fn 18 hazard cumulative gamma
*/

/*
Likelihood function for presence detector type
'simple' version is equivalent to Royle & Nichols 2003 Poisson
'integrated' alows for non-step detection function
MGE 2011-09-30
*/

#include "secr.h"

/*--------------------------------------------------------------------*/
/* define functions with third parameter z */

int par3 (int fn) {
    if ((fn==1) || (fn==3) || (fn == 5)  || (fn == 6)  || (fn == 7) || 
	(fn == 8) || (fn==10) || (fn == 11)  || (fn == 12)  || (fn == 13) || 
        (fn == 15) || (fn==17) || (fn == 18))
	return(1);
    else
	return(0);
}

/*--------------------------------------------------------------------*/
double pfn (
    int fn,
    double d2val,
    double g0,
    double sigma,
    double z,
    double miscparm [],
    double w2)
{
    double p = -1;
    fnptr hfn;
    double tmp[4];

    if (d2val > w2) 
        p = 0;
    else {
        hfn = gethfn (fn);
        tmp[0] = g0;
        tmp[1] = sigma;
        tmp[2] = z;
        tmp[3] = miscparm[0];
        p = hfn (tmp, sqrt(d2val));
    }
    return (p);
}

fnptr gethfn (int fn) 
{
    if (fn == 0)
        return(hn);
    else if (fn == 1)
        return(hr);
    else if (fn == 2)
        return(he);
    else if (fn == 3)
        return(hnc);
    else if (fn == 4)
        return(un);
    else if (fn == 5)
        return(hf);
    else if (fn == 6)
        return(hann);
    else if (fn == 7)
        return(hcln);
    else if (fn == 8)
        return(hcg);
    else if (fn == 9)
        return(hsigbin);
    else if (fn == 10)
        return(hsig);
    else if (fn == 11)
        return(hsigsph);
    else if (fn == 12)
        return(hsigsph);
    else if (fn == 14)
        return(hhn);
    else if (fn == 15)
        return(hhr);
    else if (fn == 16)
        return(hex);
    else if (fn == 17)
        return(han);
    else if (fn == 18)
        return(hcumg);
    else (error("unknown or invalid detection function"));
    return(hn);
}
/*--------------------------------------------------------------------*/

gfnptr getgfn (int fn) 
{
    if (fn == 0)
        return(ghn);
    else if (fn == 1)
        return(ghr);
    else if (fn == 2)
        return(ghe);
    else if (fn == 3)
        return(ghnc);
    else if (fn == 4)
        return(gun); /* error("uniform detection function not allowed here"); */
    else if (fn == 5)
        return(ghf);
    else if (fn == 6)
        return(gan);
    else if (fn == 7)
        return(gcln);
    else if (fn == 8)
        return(gcg);
    else if (fn == 9)
        return(gsigbin);
    else if (fn == 10)
        return(gsig);
    else if (fn == 11)
        return(gsigsph);
    else if (fn == 12)
        return(gsigSN);
    else if (fn == 13)
        return(gsigsphSN);
    else if (fn == 14)
        return(lhn);
    else if (fn == 15)
        return(lhr);
    else if (fn == 16)
        return(lex);
    else if (fn == 17)
        return(lan);
    else if (fn == 18)
        return(lcg);
    else (error("unknown or invalid detection function"));
    return(ghn);
}

gfnLptr getgfnL (int fn) 
{
    if (fn == 0)
        return(ghnL);
    else if (fn == 1)
        return(ghrL);
    else if (fn == 2)
        return(gheL);
    else if (fn == 3)
        return(ghncL);
    else if (fn == 4)
        return(gunL); /* error("uniform detection function not allowed here"); */
    else if (fn == 5)
        return(ghfL);
    else if (fn == 6)
        return(ganL);
    else if (fn == 7)
        return(gclnL);
    else if (fn == 8)
        return(gcgL);
    else if (fn == 9)
        return(gsigbinL);
    else if (fn == 10)
        return(gsigL);
    else if (fn == 11)
        return(gsigsphL);
    else if (fn == 12)
        return(gsigSNL);
    else if (fn == 13)
        return(gsigsphSNL);
    else if (fn == 14)
        return(lhnL);
    else if (fn == 15)
        return(lhrL);
    else if (fn == 16)
        return(lexL);
    else if (fn == 17)
        return(lanL);
    else if (fn == 18)
        return(lcgL);
    else (error("unknown or invalid detection function"));
    return(ghnL);
}

/*--------------------------------------------------------------------*/
double hn (double param [], double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
/*--------------------------------------------------------------------*/
double hr (double param [], double r) {
    return(param[0] * (1 - exp(- pow(r / param[1], -param[2]))));
}
/*--------------------------------------------------------------------*/
double he (double param [], double r) {
    return (param[0] * exp(-r / param[1]));
}
/*--------------------------------------------------------------------*/
double hnc (double param [], double r) {
    double temp;
    temp = param[0] * exp(- r * r  / 2 / param[1] / param[1]);
    if (round(param[2]) > 1) temp = 1 - pow(1 - temp, param[2]);
    return (temp);
}
/*--------------------------------------------------------------------*/
double un (double param [], double r) {
    if (r<param[1]) return (param[0]);
    else return (0);
}
/*--------------------------------------------------------------------*/
double hf (double param [], double r) {
    if (r<param[2]) return (param[0]);
    else return (param[0] * exp(-(r-param[2]) / param[1]));
}
/*--------------------------------------------------------------------*/
double hann (double param [], double r) {
    return (param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
        param[1] / param[1]));
}
/*--------------------------------------------------------------------*/
double hcln (double param [], double r) {
    double g0, sigma, z, CV2, meanlog, sdlog;
    g0 = param[0];
    sigma = param[1];
    z = param[2];
    CV2 = z*z/sigma/sigma;
    meanlog = log(sigma) - log(1 + CV2)/2;
    sdlog = sqrt(log(1 + CV2));
    return g0 * plnorm(r,meanlog,sdlog,0,0); 
}
/*--------------------------------------------------------------------*/
double hcg (double param [], double r) {
    return param[0] * pgamma(r,param[2],param[1]/param[2],0,0); 
}
/*--------------------------------------------------------------------*/

/* binary signal strength - (beta0-c)/sdS, beta1/sdS */
double hsigbin  (double param [], double r) {
    double gam, b0, b1;
    b0 = param[0];
    b1 = param[1];
    gam = -(b0 + b1 * r);
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength - beta0, beta1, sdS */
double hsig (double param [], double r) {
    double mu, gam;
    double beta0, beta1, sdS, cut;
    beta0 = param[0];
    beta1 = param[1];
    sdS = param[2];
    cut = param[3];
    mu = beta0 + beta1 * r;
    gam = (cut - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength with spherical spreading - beta0, beta1, sdS */
double hsigsph (double param [], double r) {
    double mu, gam;
    double beta0, beta1, sdS, cut;
    beta0 = param[0];
    beta1 = param[1];
    sdS = param[2];
    cut = param[3];
    mu =  beta0 + beta1 * (r-1) - 10 * log(r*r) / M_LN10;
    gam = (cut - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* hazard halfnormal */
double hhn (double param [], double r) {
    /* Rprintf("%6.3f %6.3f %6.3f \n", r, param[0], param[1]); */
    return(1 - exp( - param[0] * exp(- r * r / 2 / param[1] / param[1])));
}
/*--------------------------------------------------------------------*/

/* hazard hazard rate */
double hhr (double param [], double r) {
    return(1 - exp( - param[0] * ( 1 - exp(- pow(r / param[1], -param[2])))));
}
/*--------------------------------------------------------------------*/

/* hazard exponential */
double hex (double param [], double r) {
    return (1 - exp( - param[0] * exp(-r / param[1])));
}
/*--------------------------------------------------------------------*/

/* hazard annular normal */
double han (double param [], double r) {
    return (1 - exp( - param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
				   param[1] / param[1])));
}
/*--------------------------------------------------------------------*/

/* hazard cumulative gamma */
double hcumg (double param [], double r) {
    return (1 - exp( - (1 - exp( - param[0] * exp(-r / param[1])))));
}
/*====================================================================*/

/* halfnormal */
double ghn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    return (gsbval[c] * exp(-d2(k, m, traps, mask, kk, mm) / 2 /
        gsbval[cc + c] / gsbval[cc + c]));
}
/*--------------------------------------------------------------------*/

/* compound halfnormal */
double ghnc
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double temp;
    temp = exp(-d2(k, m, traps, mask, kk, mm) / 2 / gsbval[cc + c] / gsbval[cc + c]);
    temp = 1 - pow(1 - temp, gsbval[2*cc + c]);

    return (gsbval[c] * temp);
}
/*--------------------------------------------------------------------*/

/* hazard rate */
double ghr
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    return (gsbval[c] * (1 - exp(-
        pow(sqrt(d2(k,m,traps,mask,kk,mm)) / gsbval[cc + c], - gsbval[cc * 2 + c]))));
}
/*--------------------------------------------------------------------*/

/* exponential */
double ghe
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    return (gsbval[c] * exp(-sqrt(d2(k,m,traps,mask,kk,mm)) / gsbval[cc + c]));
}
/*--------------------------------------------------------------------*/

/* uniform */
double gun
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    if (d2(k, m, traps, mask, kk, mm) < (gsbval[cc + c]*gsbval[cc + c]))
	return (gsbval[c]);
    else return (0);
}
/*--------------------------------------------------------------------*/

/* 'flat-topped exponential' 2009 09 01 */
double ghf
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, w, g0, sigma;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    w = gsbval[cc * 2 + c];
    if (d<w) return (g0);
    else return (g0 * exp(-(d-w) / sigma));
}
/*--------------------------------------------------------------------*/

/* annular halfnormal 2010-06-15 */
double gan
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, w, g0, sigma;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    w = gsbval[cc * 2 + c];
    return (g0 * exp(-(d-w)*(d-w) / 2 / sigma / sigma));
}
/*--------------------------------------------------------------------*/

/* cumulative lognormal 2010-10-10 */
double gcln
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, g0, sigma, z, CV2, meanlog, sdlog;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    CV2 = z*z/sigma/sigma;
    meanlog = log(sigma) - log(1 + CV2)/2;
    sdlog = sqrt(log(1 + CV2));
    return g0 * plnorm(d,meanlog,sdlog,0,0); 
}
/*--------------------------------------------------------------------*/

/* cumulative normal 2010-10-11 */
double gcn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, g0, sigma, z;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return g0 * pnorm(d,sigma,z,0,0); 
}
/*--------------------------------------------------------------------*/

/* cumulative gamma 2010-10-13 */
double gcg
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, g0, sigma, z;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return g0 * pgamma(d,z,sigma/z,0,0); 
}
/*--------------------------------------------------------------------*/

/* reverse sigmoid */
double grs
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, g0, sigma, z, x;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    x = z * (d - sigma);
    return( g0 * (1 + (1 - expmin(x)) / (1 + expmin(x))) / 2 ) ;
}
/*--------------------------------------------------------------------*/

/* binary signal strength - (beta0-c)/sdS, beta1/sdS */
double gsigbin
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double gam, b0, b1;
    b0 = gsbval[c];
    b1 = gsbval[cc + c];
    gam = -(b0 + b1 * sqrt(d2(k,m, traps, mask, kk, mm)));
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength - beta0, beta1, sdS */
double gsig
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm []) {
    double mu, gam, sdS;
    mu = mufn (k, m, gsbval[c], gsbval[cc + c], traps, mask, kk, mm, 0);
    sdS = gsbval[cc * 2 + c];
    gam = (miscparm[0] - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength with spherical spreading - beta0, beta1, sdS */
double gsigsph
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm []) {
    double mu, gam, sdS;
    mu = mufn (k, m, gsbval[c], gsbval[cc + c], traps, mask, kk, mm, 1);
    sdS = gsbval[cc * 2 + c];
    gam = (miscparm[0] - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal-noise without spherical spreading - beta0, beta1, sdS */
double gsigSN
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm []) {
    double cut, muS, sdS, muN, sdN;
    muS = mufn (k, m, gsbval[c], gsbval[cc + c], traps, mask, kk, mm, 0);
    sdS = gsbval[cc * 2 + c];
    cut = miscparm[0];
    muN = miscparm[1];
    sdN = miscparm[2];
    return (pnorm(cut, muS-muN, sqrt(sdS*sdS+sdN*sdN), 0, 0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal-noise with spherical spreading - beta0, beta1, sdS */
/* interpret sdS as sqrt(var(S-N)) */
double gsigsphSN
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm []) {
    double cut, muS, sdS, muN, sdN;
    muS = mufn (k, m, gsbval[c], gsbval[cc + c], traps, mask, kk, mm, 1);
    sdS = gsbval[cc * 2 + c];
    cut = miscparm[0];
    muN = miscparm[1];    
    sdN = miscparm[2];    
    return (pnorm(cut, muS-muN, sqrt(sdS*sdS+sdN*sdN), 0, 0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* hazard halfnormal */
double lhn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double lambda0, sigma;
    double invdenom = 1.0;
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
/*
    if (fabs(miscparm[0]) > 0.5) invdenom = mask[2*mm+m];
*/
    return (1 - exp(- lambda0 * invdenom *  exp(-d2(k, m, traps, mask, kk, mm)
						/ 2 / sigma / sigma)));    
}
/*--------------------------------------------------------------------*/

/* hazard hazard rate */
double lhr
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, lambda0, sigma, z;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return (1 - exp(- lambda0 * ( 1 - exp(- pow(d /sigma , - z)))));
}
/*--------------------------------------------------------------------*/

/* hazard exponential */
double lex
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double lambda0, sigma;
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    return (1 - exp(-lambda0 * exp(-sqrt(d2(k,m,traps,mask,kk,mm))
						/ sigma)));
}
/*--------------------------------------------------------------------*/

/* hazard annular normal */
double lan
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, w, lambda0, sigma;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    w = gsbval[cc * 2 + c];
    return (1 - exp(-lambda0 * exp(-(d-w)*(d-w) / 2 / sigma / sigma)));
}
/*--------------------------------------------------------------------*/

/* hazard cumulative gamma */
double lcg
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm [])
{
    double d, lambda0, sigma, z;
    d = sqrt(d2(k,m,traps,mask,kk,mm));
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return (1 - exp(- lambda0 * pgamma(d,z,sigma/z,0,0))); 
}
/*====================================================================*/

/* halfnormal */
double ghnL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    return (gsbval[c] * exp(-d2L(k, m, dist2, kk) / 2 /
        gsbval[cc + c] / gsbval[cc + c]));
}
/*--------------------------------------------------------------------*/

/* compound halfnormal */
double ghncL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double temp;
    temp = exp(-d2L(k, m, dist2, kk) / 2 / gsbval[cc + c] / gsbval[cc + c]);
    temp = 1 - pow(1 - temp, gsbval[2*cc + c]);

    return (gsbval[c] * temp);
}
/*--------------------------------------------------------------------*/

/* hazard rate */
double ghrL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    return (gsbval[c] * (1 - exp(-
        pow(sqrt(d2L(k, m, dist2, kk)) / gsbval[cc + c], - gsbval[cc * 2 + c]))));
}
/*--------------------------------------------------------------------*/

/* exponential */
double gheL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    return (gsbval[c] * exp(-sqrt(d2L(k, m, dist2, kk)) / gsbval[cc + c]));
}
/*--------------------------------------------------------------------*/

/* uniform */
double gunL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    if (d2L(k, m, dist2, kk) < (gsbval[cc + c]*gsbval[cc + c]))
	return (gsbval[c]);
    else return (0);
}
/*--------------------------------------------------------------------*/

/* 'flat-topped exponential' 2009 09 01 */
double ghfL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, w, g0, sigma;
    d = sqrt(d2L(k, m, dist2, kk));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    w = gsbval[cc * 2 + c];
    if (d<w) return (g0);
    else return (g0 * exp(-(d-w) / sigma));
}
/*--------------------------------------------------------------------*/

/* annular halfnormal 2010-06-15 */
double ganL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, w, g0, sigma;
    d = sqrt(d2L(k, m, dist2, kk));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    w = gsbval[cc * 2 + c];
    return (g0 * exp(-(d-w)*(d-w) / 2 / sigma / sigma));
}
/*--------------------------------------------------------------------*/

/* cumulative lognormal 2010-10-10 */
double gclnL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, g0, sigma, z, CV2, meanlog, sdlog;
    d = sqrt(d2L(k, m, dist2, kk));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    CV2 = z*z/sigma/sigma;
    meanlog = log(sigma) - log(1 + CV2)/2;
    sdlog = sqrt(log(1 + CV2));
    return g0 * plnorm(d,meanlog,sdlog,0,0); 
}
/*--------------------------------------------------------------------*/

/* cumulative normal 2010-10-11 */
double gcnL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, g0, sigma, z;
    d = sqrt(d2L(k, m, dist2, kk));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return g0 * pnorm(d,sigma,z,0,0); 
}
/*--------------------------------------------------------------------*/

/* cumulative gamma 2010-10-13 */
double gcgL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, g0, sigma, z;
    d = sqrt(d2L(k, m, dist2, kk));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return g0 * pgamma(d,z,sigma/z,0,0); 
}
/*--------------------------------------------------------------------*/

/* reverse sigmoid */
double grsL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, g0, sigma, z, x;
    d = sqrt(d2L(k, m, dist2, kk));
    g0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    x = z * (d - sigma);
    return( g0 * (1 + (1 - expmin(x)) / (1 + expmin(x))) / 2 ) ;
}
/*--------------------------------------------------------------------*/

/* binary signal strength - (beta0-c)/sdS, beta1/sdS */
double gsigbinL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double gam, b0, b1;
    b0 = gsbval[c];
    b1 = gsbval[cc + c];
    gam = -(b0 + b1 * sqrt(d2L(k, m, dist2, kk)));
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength - beta0, beta1, sdS */
double gsigL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm []) {
    double mu, gam, sdS;
    mu = mufnL (k, m, gsbval[c], gsbval[cc + c], dist2, kk, 0);
    sdS = gsbval[cc * 2 + c];
    gam = (miscparm[0] - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength with spherical spreading - beta0, beta1, sdS */
double gsigsphL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm []) {
    double mu, gam, sdS;
    mu = mufnL (k, m, gsbval[c], gsbval[cc + c], dist2, kk, 1);
    sdS = gsbval[cc * 2 + c];
    gam = (miscparm[0] - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal-noise without spherical spreading - beta0, beta1, sdS */
double gsigSNL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm []) {
    double cut, muS, sdS, muN, sdN;
    muS = mufnL (k, m, gsbval[c], gsbval[cc + c], dist2, kk, 0);
    sdS = gsbval[cc * 2 + c];
    cut = miscparm[0];
    muN = miscparm[1];
    sdN = miscparm[2];
    return (pnorm(cut, muS-muN, sqrt(sdS*sdS+sdN*sdN), 0, 0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal-noise with spherical spreading - beta0, beta1, sdS */
/* interpret sdS as sqrt(var(S-N)) */
double gsigsphSNL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm []) {
    double cut, muS, sdS, muN, sdN;
    muS = mufnL (k, m, gsbval[c], gsbval[cc + c], dist2, kk, 1);
    sdS = gsbval[cc * 2 + c];
    cut = miscparm[0];
    muN = miscparm[1];    
    sdN = miscparm[2];    
    return (pnorm(cut, muS-muN, sqrt(sdS*sdS+sdN*sdN), 0, 0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* hazard halfnormal */
double lhnL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double lambda0, sigma;
    double invdenom = 1.0;
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    /* if (fabs(miscparm[0]) > 0.5) invdenom = mask[2*mm+m]; */
    return (1 - exp(- lambda0 * invdenom *  exp(-d2L(k, m, dist2, kk)
						/ 2 / sigma / sigma)));    
}
/*--------------------------------------------------------------------*/

/* hazard hazard rate */
double lhrL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, lambda0, sigma, z;
    double invdenom = 1.0;
    /* if (fabs(miscparm[0]) > 0.5) invdenom = mask[2*mm+m]; */
    d = sqrt(d2L(k, m, dist2, kk));
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return (1 - exp(- lambda0 * invdenom * ( 1 - exp(- pow(d /sigma , - z)))));
}
/*--------------------------------------------------------------------*/

/* hazard exponential */
double lexL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double lambda0, sigma;
    double invdenom = 1.0;
    /* if (fabs(miscparm[0]) > 0.5) invdenom = mask[2*mm+m]; */
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    return (1 - exp(-lambda0 * invdenom * exp(-sqrt(d2L(k, m, dist2, kk))
						/ sigma)));
}
/*--------------------------------------------------------------------*/

/* hazard annular normal */
double lanL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, w, lambda0, sigma;
    double invdenom = 1.0;
    /* if (fabs(miscparm[0]) > 0.5) invdenom = mask[2*mm+m]; */
    d = sqrt(d2L(k, m, dist2, kk));
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    w = gsbval[cc * 2 + c];
    return (1 - exp(-lambda0 * invdenom * exp(-(d-w)*(d-w) / 2 / sigma / sigma)));
}
/*--------------------------------------------------------------------*/

/* hazard cumulative gamma */
double lcgL
    (int k, int m, int c, double gsbval[], int cc, double dist2[],
    int kk, double miscparm [])
{
    double d, lambda0, sigma, z;
    double invdenom = 1.0;
    /* if (fabs(miscparm[0]) > 0.5) invdenom = mask[2*mm+m]; */
    d = sqrt(d2L(k, m, dist2, kk));
    lambda0 = gsbval[c];
    sigma = gsbval[cc + c];
    z = gsbval[cc * 2 + c];
    return (1 - exp(- lambda0 * invdenom * pgamma(d,z,sigma/z,0,0))); 
}
/*--------------------------------------------------------------------*/

double mufn (
    int k,
    int m,
    double b0,
    double b1,
    double A1[],
    double A2[],
    int A1rows,
    int A2rows,
    int spherical)
/*
   Return predicted signal strength at m for source at point k,
   given strength at source of b0 dB and attenuation of b1 dB/m.
   Spherical spreading is included if spherical > 0
   Coordinates of points are in A1 and A2 which have respectively
   A1rows and A2rows
*/
{
    double d2val;
    d2val = d2(k,m, A1, A2, A1rows, A2rows);
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
/*==============================================================================*/

double mufnL (
    int k,
    int m,
    double b0,
    double b1,
    double dist2[],
    int kk,
    int spherical)
/*
   Return predicted signal strength at k for source at point m,
   given strength at source of b0 dB and attenuation of b1 dB/m.
   Spherical spreading is included if spherical > 0
   Uses distance lookup in dist2
*/
{
    double d2val;
    d2val = d2L(k, m, dist2, kk);
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
/*==============================================================================*/

