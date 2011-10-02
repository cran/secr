/*
    Detection functions
    2011-09-30
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
*/

/*
Likelihood function for presence detector type
'simple' version is equivalent to Royle & Nichols 2003 Poisson
'integrated' alows for non-step detection function
MGE 2011-09-30
*/

#include "secr.h"

double pfn (
    int fn,
    double d2val,
    double g0,
    double sigma,
    double z,
    double cut,
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
        tmp[3] = cut;
        p = hfn (tmp, sqrt(d2val));
    }
    return (p);
}

fnptr gethfn (int fn) 
{
    if (fn == 0)
        return(hn);
    else if (fn == 1)
        return(hz);
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
    else (error("unknown or invalid detection function"));
    return(hn);
}
/*--------------------------------------------------------------------*/

gfnptr getgfn (int fn) 
{
    if (fn == 0)
        return(ghn);
    else if (fn == 1)
        return(ghz);
    else if (fn == 2)
        return(ghe);
    else if (fn == 3)
        return(ghnc);
    else if (fn == 4)
        error("uniform detection function not allowed here");
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
    else (error("unknown or invalid detection function"));
    return(ghn);
}

/*--------------------------------------------------------------------*/
double hn (double param [], double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
/*--------------------------------------------------------------------*/
double hz (double param [], double r) {
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

/* halfnormal */
double ghn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double cut)
{
    return (gsbval[c] * exp(-d2(k, m, traps, mask, kk, mm) / 2 /
        gsbval[cc + c] / gsbval[cc + c]));
}
/*--------------------------------------------------------------------*/

/* compound halfnormal */
double ghnc
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double cut)
{
    double temp;
    temp = exp(-d2(k, m, traps, mask, kk, mm) / 2 / gsbval[cc + c] / gsbval[cc + c]);
    temp = 1 - pow(1 - temp, gsbval[2*cc + c]);

    return (gsbval[c] * temp);
}
/*--------------------------------------------------------------------*/

/* hazard rate */
double ghz
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double cut)
{
    return (gsbval[c] * (1 - exp(-
        pow(sqrt(d2(k,m,traps,mask,kk,mm)) / gsbval[cc + c], - gsbval[cc * 2 + c]))));
}
/*--------------------------------------------------------------------*/

/* exponential */
double ghe
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double cut)
{
    return (gsbval[c] * exp(-sqrt(d2(k,m,traps,mask,kk,mm)) / gsbval[cc + c]));
}
/*--------------------------------------------------------------------*/

/* 'flat-topped exponential' 2009 09 01 */
double ghf
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut)
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
    double mask[], int kk, int mm, double cut) {
    double mu, gam, sdS;
    mu = mufn (k, m, gsbval[c], gsbval[cc + c], traps, mask, kk, mm);
    sdS = gsbval[cc * 2 + c];
    gam = (cut - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

/* signal strength with spherical spreading - beta0, beta1, sdS */
double gsigsph
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double cut) {
    double mu, gam, sdS;
    mu = mufnsph (k, m, gsbval[c], gsbval[cc + c],
         traps, mask, kk, mm);
    sdS = gsbval[cc * 2 + c];
    gam = (cut - mu) / sdS;
    return (pnorm(gam,0,1,0,0));    /* upper */
}
/*--------------------------------------------------------------------*/

