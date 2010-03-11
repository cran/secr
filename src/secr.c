/* 
   External procedures for secr package
   Murray Efford 2008 07 10 - 2009 12 12
   2009 09 10 single-catch trap simulation integrated with simsecr!
   2009 09 24 remove local math functions & rely on R
   2009 09 26 convert all comments to C (not C++ double slash)
   2009 09 26 un-mix declarations and code in integralprw1 & secrloglik
   2009 09 30 remaining code incompatible with ISO C: variable-length arrays  
   2009 09 30 amendments by Ray Brownrigg to compile on Unix
              - removed all mention of EXPORT
              - changed random() to Random  
   ...
   2009 10 25 MR secr code tweaked
   2009 10 26 Major revision to include count, area and signal detectors
   2009 11 .. Polygon detector added
   2009 11 21 allocate 'capt' on stack with R_alloc in trapping simulation routines
   2009 12 03 completed 'transect' detector and revision of 'signal' detector
   2009 12 04 combined proximity, count, quadrat detectors
   2009 12 04 rescale parameter of binomial: detection function now models Np i.e. E(x), just like Poisson & Bernoulli
   2009 12 12 finite mixtures completed
   2010 03 09 fixed bug in rdiscrete
 
   can compile with gcc 4.2.1-sjlj :
   gcc -Ic:/progra~1/R/R-2.9.2/include -c secr.c -Wall -pedantic -std=gnu99
*/

#include "secr.h"
#include <math.h>    
#include <stdlib.h>
#include <stdio.h> 
#include <R.h>       /* random numbers */
#include <Rmath.h>   /* R math functions e.g. dbinom, dpois */
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/*==============================================================================*/
/* function pointers */

typedef double (*gfnptr)(int, int, int, double[], int, double[], double[], int, int, double);
/*
    (int k, int m, int c, int double gsbval[], int cc, double traps[], 
     double mask[], int kk, int mm, double cut);
*/

typedef double (*prwfnptr)(int, int, int, int, int, int[], double[], double[], int[], 
    double[], int, double[], int, int, int, int, int, int, gfnptr, double[], double[], 
    double[], double);
/*
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[], 
     double gk[], int binomN, double detspec[], int cc, int nc, int kk, int ss,
     int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], double minp);
*/

typedef double (*fnptr)(double[], double);

/*==============================================================================*/

double huge = 1e10;
double fuzz = 1e-200;
double tiny = 1e-6;         /* added to ensure reliable typecast of double to integer */
double minimumexp = -100;
int maxnpoly = 1000;       /* applies to polygons & transects */
int maxperpoly = 100;
int maxnmix = 2;          /* maximum number of mixtures */

struct trap_animal {                   
    int     trap;
    int     animal;
    double  time;
};
struct rpoint {                   
    double x;
    double y;
};

FILE *out;      /* for debugging */
/*==============================================================================*/

void R_CheckUserInterrupt(void);

/*==============================================================================*/

double Random () { 
/*
   this call to R requires preceding 
   GetRNGstate();
   and following
   PutRNGstate();
*/
    return( unif_rand() );     /* 2009 09 08 */
}
/*==============================================================================*/

double expmin (double x)
{
  if (x < minimumexp)
      return(0);
  else 
      return(exp(x));  
}
/*==============================================================================*/

/* index to vector element corresponding to cell i,j,k in 3D array
   stored in column-major order */

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}
/*==============================================================================*/

/* index to vector element corresponding to cell i,j,k,l in 4D array
   stored in column-major order */

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}
/*==============================================================================*/

double gpois (int count, double lambda, int uselog)
{  
    if (count == 0) {
        if (uselog) 
            return (-lambda);
        else 
            return (exp(-lambda));  
    }
    else
        return (dpois(count, lambda, uselog));        
}

double gbinom (int count, int size, double p, int uselog)
{  
    double x;
    int i; 
    if (count == 0) {
        p = 1 - p;
        x = p;
        for (i=1; i< size; i++) x = x*p;
        if (uselog) x = log(x);
        return (x);   /* faster */
    } 
    else return (dbinom (count, size, p, uselog));
}

double gnbinom (int count, int size, double mu, int uselog)
{  
    /* prob = size / (size + mu) */
    size = abs(size);  /* in case negative 'binomN' passed */
    return (dnbinom (count, size, size/(size+mu), uselog));
}

/*==============================================================================*/

double d2 (
    int k, 
    int m, 
    double A1[], 
    double A2[], 
    int A1rows, 
    int A2rows)
/*
   return squared distance between two points given by row k in A1
   and row m in A2, where A1 and A2 have respectively A1rows and A2rows
*/
{
    return(
        (A1[k] - A2[m]) * (A1[k] - A2[m]) +
        (A1[k + A1rows] - A2[m + A2rows]) * (A1[k + A1rows] - A2[m + A2rows])
    );
}

/*==============================================================================*/

double mufn (
    int k, 
    int m, 
    double b0, 
    double b1,
    double A1[], 
    double A2[], 
    int A1rows, 
    int A2rows)
/*
   Return predicted signal strength at m for source at point k, 
   given strength at source of b0 dB and attenuation of b1 dB/m.
   Spherical spreading is NOT included
   Coordinates of points are in A1 and A2 which have respectively 
   A1rows and A2rows
*/
{
    double d2val;
    d2val = d2(k,m, A1, A2, A1rows, A2rows);  
    return (b0 + b1 * sqrt(d2val));
}
/*==============================================================================*/

double mufnsph (
    int k, 
    int m, 
    double b0, 
    double b1,
    double A1[], 
    double A2[], 
    int A1rows, 
    int A2rows)
/*
   Return predicted signal strength at m for source at point k, 
   given strength at source of b0 dB and attenuation of b1 dB/m.
   Spherical spreading is included
   Coordinates of points are in A1 and A2 which have respectively 
   A1rows and A2rows
*/
{
    double d2val;
    d2val = d2(k,m, A1, A2, A1rows, A2rows);  
     if (d2val>1) {            
        return (b0 - 10 * log ( d2val ) / 2.302585 + b1 * (sqrt(d2val)-1));  /* checked 2008 08 08 */
    }
    else 
        return (b0);
}
/*==============================================================================*/

double hn (double param [], double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));  
}
double hz (double param [], double r) {
    return(param[0] * (1 - exp(- pow(r / param[1], -param[2]))));  
}
double he (double param [], double r) { 
    return (param[0] * exp(-r / param[1]));  
}
double hnc (double param [], double r) {
    double temp;
    temp = param[0] * exp(- r * r  / 2 / param[1] / param[1]);  
    if (round(param[2]) > 1) temp = 1 - pow(1 - temp, param[2]);
    return (temp);
}
double un (double param [], double r) { 
    if (r<param[1]) return (param[0]);
    else return (0);  
}
double hf (double param [], double r) { 
    if (r<param[2]) return (param[0]);
    else return (param[0] * exp(-(r-param[2]) / param[1]));  
}

/* halfnormal */
double ghn 
    (int k, int m, int c, double gsbval[], int cc, double traps[], 
    double mask[], int kk, int mm, double cut)
{
    return (gsbval[c] * exp(-d2(k, m, traps, mask, kk, mm) / 2 / 
        gsbval[cc + c] / gsbval[cc + c]));  
}

/* compound halfnormal */
double ghnc 
    (int k, int m, int c, double gsbval[], int cc, double traps[], 
    double mask[], int kk, int mm, double cut)
{
    double temp;
    temp = gsbval[c] * exp(-d2(k, m, traps, mask, kk, mm) / 2 / gsbval[cc + c] / 
        gsbval[cc + c]);  
    if (gsbval[2*cc + c] > 1) temp = 1 - pow(1 - temp, gsbval[2*cc + c]); 
    return (temp);
}

/* hazard rate */
double ghz 
    (int k, int m, int c, double gsbval[], int cc, double traps[], 
    double mask[], int kk, int mm, double cut) 
{
    return (gsbval[c] * (1 - exp(- 
        pow(sqrt(d2(k,m,traps,mask,kk,mm)) / gsbval[cc + c], - gsbval[cc * 2 + c])))); 
}
/* exponential */
double ghe 
    (int k, int m, int c, double gsbval[], int cc, double traps[], 
    double mask[], int kk, int mm, double cut) 
{
    return (gsbval[c] * exp(-sqrt(d2(k,m,traps,mask,kk,mm)) / gsbval[cc + c]));  
}
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

/*===============================================================*/
void rgr(double *x, int n, void *ex) {
    int i;
    int fn;
    double * p;
    double tmp[4];
    p = (double*) ex;
    for (i=0; i<4; i++) tmp[i] = p[i];
    fn = tmp[3];
    fnptr fnp; 
    fnp = hn;    /* default */
    if (fn == 1)
        fnp = hz;
    else if (fn == 2)
        fnp = he;
    else if (fn == 3)
        fnp = hnc;
    else if (fn == 4)
        fnp = un;
    else if (fn == 5)
        fnp = hf;   
    for (i=0; i<n; i++) {
        x[i] = x[i] * fnp(tmp,x[i]);   /* r.g(r) */
    }
}

void gxy (int *n, int *fn, double *par, double *w, double *xy) {
    int maxj = 100000;
    double r;
    double theta;
    fnptr fnp; 
    int i = 0;
    int j;
    fnp = hn;
    if (*fn == 1)
        fnp = hz;
    else if (*fn == 2)
        fnp = he;
    else if (*fn == 3)
        fnp = hnc;
    else if (*fn == 4)
        fnp = un;
    else if (*fn == 5)
        fnp = hf;
    for (i=0; i< *n; i++) {
        for (j=0; j<maxj; j++) {
            r = *w * sqrt(unif_rand());
            if (unif_rand() < fnp(par, r))
                break;
        }
        theta = unif_rand() * 2 * M_PI;
        xy[i]      = r * cos(theta);
        xy[*n + i] = r * sin(theta);
    }
}

double distance (struct rpoint p1, struct rpoint p2) {
    return(sqrt ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y)));
}
double gr (int *fn, double par[], struct rpoint xy, struct rpoint animal) {
    double r;
    fnptr fnp; 
    fnp = hn;
    if (*fn == 1)
        fnp = hz;
    else if (*fn == 2)
        fnp = he;
    else if (*fn == 3)
        fnp = hnc;
    else if (*fn == 4)
        fnp = un;
    else if (*fn == 5)
        fnp = hf;
    r = distance (xy, animal);
    return (fnp(par,r));
}

double gintegral (int fn, double par[]) {
/* integral of radial 2-D function */
    double ex[4];  
    double a;
    int b;
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
    a = 0;
    b = 1;    /* signals bounds 0,Inf */
    ex[0] = par[0];
    ex[1] = par[1];
    ex[2] = par[2];
    ex[3] = fn;

    Rdqagi(rgr, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);

    /* ignoring ier etc. */
    return (result * 2 * M_PI);
}

/* find upper and lower points on perimeter of poly at x-coordinate x */
void yab(double x[], int *i, int *np, double poly[], double *a, double *b) {
    int k;
    int nv = 0;
    double ab[3];
    /* note 'sign' is RMath function */
    for (k=0; k< (*np-1); k++) {
        if (sign(poly[k]- x[*i]) != sign(poly[k+1]- x[*i])) {
           ab[nv] = poly[k+ *np] + (x[*i]-poly[k]) * (poly[k+1+*np]-poly[k+*np]) / (poly[k+1]-poly[k]);
           nv++;
        }
        if (nv>2) break;
    }
    if (ab[0]>ab[1])
        {*a = ab[1]; *b = ab[0];}
    else
        {*a = ab[0]; *b = ab[1];}

}

void fy(double *x, int n, void *ex) {
    int i;
    int fn;
    double * p;
    double mx,my;
    double xy[2];
    double d;
    p = (double*) ex;
    fn = round(p[3]);
    mx = p[4];
    my = p[5];
    xy[0] = p[6];

    /* set detection function */
    fnptr fnp; 
    fnp = hn;    /* default */
    if (fn == 1)
        fnp = hz;
    else if (fn == 2)
        fnp = he;
    else if (fn == 3)
        fnp = hnc;
    else if (fn == 4)
        fnp = un;
    else if (fn == 5)
        fnp = hf;   
    for (i=0; i<n; i++) {
        xy[1] = x[i];   /* set each y value */
        d = sqrt ( (xy[1]-my)*(xy[1]-my) + (xy[0]-mx)*(xy[0]-mx) );
        x[i] = fnp(p, d);   /* g(r) */
    }
}

void fx(double *x, int n, void *ex) {
    int i;
    double * p;
    double * a;
    double * b;
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
    double *poly;
    int kk;
    p = (double*) ex;
    kk = round(p[9]);
    a = (double *) R_alloc(1, sizeof(double));
    b = (double *) R_alloc(1, sizeof(double));
    poly = (double *) R_alloc(kk * 2, sizeof(double));
    for (i=0; i<kk; i++) {
        poly[i] = p[i+10];
        poly[i+kk] = p[i+kk+10];
    }
    for (i=0; i<n; i++) {
        yab(x, &i, &kk, poly, a, b);   /* refine limits here */
        p[6] = x[i];                   /* pass forward the value of x; consider &ex etc. */
        Rdqags(fy, ex, a, b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
        x[i] = result;
    }
}

double integral2D  (int fn, int m, int c, double gsbval[], int cc, double traps[], 
    double mask[], int n1, int n2, int kk, int mm) {
    double *ex; 
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

    /* limits from bounding box of this polygon */
    ns = n2-n1+1;
    for (k=0; k<ns; k++) {
        ax = fmin2(ax, traps[k+n1]);
        bx = fmax2(bx, traps[k+n1]);
    }

    /* pass parameters etc. through pointer */
    ex = (double *) R_alloc(10 + 2 * ns, sizeof(double));
    ex[0] = gsbval[c];    
    ex[1] = gsbval[cc + c];
    ex[2] = gsbval[2*cc + c];
    ex[3] = fn;
    ex[4] = mask[m];
    ex[5] = mask[m+mm];
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;

    /* also pass polygon vertices */
    for (k=0; k<ns; k++) {
        ex[k+10] = traps[k+n1];        /* x */
        ex[k+ns+10] = traps[k+n1+kk];  /* y */
    }
    Rdqags(fx, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
    if (ier != 0) Rprintf("ier error code in integral2D %5d\n", ier);
    return (result); 
}

/*===============================================================*/
void integral2Dtest
    (int *fn, int *m, int *c, double *gsbval, int *cc, double *traps, 
    double *mask, int *n1, int *n2, int *kk, int *mm, double *result)
{
    double *ex; 
    double ax=1e20;
    double bx=-1e20;
    double res;
    double epsabs = 0.0001;
    double epsrel = 0.0001;
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
    double par[4];

    /* limits from bounding box of this polygon */
    ns = *n2 - *n1 + 1;
    for (k=0; k<ns; k++) {
        ax = fmin2(ax, traps[k+ns]);
        bx = fmax2(bx, traps[k+ns]);
    }

    par[0] = 1;
    par[1] = gsbval[*cc + *c];
    par[2] = gsbval[2* *cc + *c];   

    ex = (double *) R_alloc(10 + 2 * *kk, sizeof(double));
    ex[0] = gsbval[*c];    
    ex[1] = gsbval[*cc + *c];
    ex[2] = gsbval[2* *cc + *c];
    ex[3] = *fn;
    ex[4] = mask[*m];
    ex[5] = mask[*m+ *mm];
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;

    for (k=0; k<ns; k++) {
        ex[k+10] = traps[k+ *n1];
        ex[k+ns+10] = traps[k+ *n1 + *kk];
    }
    Rdqags(fx, ex, &ax, &bx, &epsabs, &epsrel, &res, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
    *result = res;
}

/*===============================================================*/

struct rpoint getxy(double l, double cumd[], struct rpoint line[], int kk, double offset) {
/* return the xy coordinates of point l metres along a transect */
/* offset is the starting position for this transect */
    int k;
    double pr, d, d12;
    struct rpoint xy;
    for (k=offset+1; k<(offset+kk); k++) if (cumd[k]>l) break;
    d = l - cumd[k-1];  /* distance along leg */
    d12 = cumd[k] - cumd[k-1];
    if (d12>0)
        pr = d / d12;
    else
        pr = 0;
    xy.x = line[k-1].x + (line[k].x - line[k-1].x) * pr;
    xy.y = line[k-1].y + (line[k].y - line[k-1].y) * pr;
    return(xy);
}

void fx1 (double *x, int n, void *ex) {
    int i;
    int ns;
    int fn;
    struct rpoint *line;
    struct rpoint mxy;
    struct rpoint xy;
    double * p;
    double *cumd;
    double d;
    fnptr fnp; 
    /* extract parameters passed in void pointer ex */
    p = (double*) ex;
    fn = round(p[3]);
    mxy.x = p[4];
    mxy.y = p[5];
    ns = round(p[9]);
    /* coordinates of vertices */
    line = (struct rpoint *) R_alloc(ns, sizeof(struct rpoint));
    for (i=0; i<ns; i++) {
        line[i].x = p[i+10];
        line[i].y = p[i+ns+10];
    }
    /* cumulative distance along line */
    cumd = (double *) R_alloc(ns + 1, sizeof(double));
    cumd[0] = 0;
    for (i=0; i<(ns-1); i++) {
        cumd[i+1] = cumd[i] + distance (line[i],line[i+1]);
    }
    /* set detection function - default hn */
    fnp = hn;         
    if (fn == 1)
        fnp = hz;
    else if (fn == 2)
        fnp = he;
    else if (fn == 3)
        fnp = hnc;
    else if (fn == 4)
        fnp = un;
    else if (fn == 5)
        fnp = hf;   
    /* for each x in x[] */
    for (i=0; i<n; i++) {              
        xy = getxy (x[i], cumd, line, ns, 0);
        d = distance (xy, mxy);
        x[i] = fnp(p, d);   /* g(r) */
    }
}

double integral1D
    (int fn, int m, int c, double gsbval[], int cc, double traps[], 
    double mask[], int n1, int n2, int kk, int mm)
{
    double *ex; 
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
    for (k=n1+1; k<=n2; k++) {         /* upper bound is length of this transect */
        bx += sqrt( (traps[k] - traps[k-1]) * (traps[k] - traps[k-1]) +
                    (traps[k+kk] - traps[k-1+kk]) * (traps[k+kk] - traps[k-1+kk]) );
    }
    /* pass parameters etc. through pointer */
    ex = (double *) R_alloc(10 + 2 * ns, sizeof(double));
    ex[0] = gsbval[c];    
    ex[1] = gsbval[cc + c];
    ex[2] = gsbval[2*cc + c];
    ex[3] = fn;
    ex[4] = mask[m];
    ex[5] = mask[m+mm];
    ex[6] = 0;
    ex[7] = 0;
    ex[8] = 0;
    ex[9] = ns;  
    for (k=0; k<ns; k++) {             /* pass transect vertices */
        ex[k+10] = traps[k+n1];        /* x */
        ex[k+ns+10] = traps[k+n1+kk];  /* y */
    }
    Rdqags(fx1, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
    if (ier != 0) Rprintf("ier error code in integral1D %5d\n", ier);
    return (result); 
}

/*===============================================================*/

double pndot (int m, int n, int s1, int s2, int x, int ncol, int gsb0[], 
    double gk0[], int ss, int kk, int cc0, int nmix)
/*
    probability animal at point m on mask is caught
    n may indicate group (full likelihood; ncol= number of groups) or 
    individual (conditional likelihood; ncol= number of individuals)
    aligned with secrloglik 2009 06 25

    2009 10 24 adjusted to allow summation over qq < ss 
    2009 11 12 'kk' should be number of parts for polygon detectors
*/
{
    int k,s,c;
    int gi;
    int wxi;
    double pp;
    pp = 1;
    for (k=0; k< kk; k++) {
        for (s=s1-1; s<s2; s++) { 
            wxi = i4(n,s,k,x,ncol,ss,kk);
            c = gsb0[wxi] - 1; 
            if (c >= 0) {    /* drops unset traps */
                gi = i3(c,k,m,cc0,kk);
                pp *= 1 - gk0[gi];
            }
        }
    }
    return (1 - pp);
}
/*===============================================================*/

void integralprw1 (
    int    *detect,    /* detector 0 multi, 1 proximity etc. */
    double *gsb0val,   /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *nc,        /* number of individuals */   
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    int    *nmix,      /* number of mixtures */ 
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *cc0,       /* number of g0/sigma/b combinations [naive animal] */
    int    *gsb0,      /* lookup which g0/sigma/b combination to use for given n, S, K [naive animal] */
    int    *ncol,      /* number of columns in gsb0; added 2009 06 25 */
    double *area,      /* area associated with each mask point (ha) */
    int    *fn,        /* codes
                          0 = halfnormal, 
                          1 = hazard rate, 
                          2 = exponential, 
                          9 = binary signal strength, 
                          10 = signal strength, 
                          11 = signal strength spher, */
    int    *binomN,    /* number of trials for 'count' detector modelled with binomial */
    double *cut,       /* transformed signal strength threshold for detection */
    double *a,         /* return value integral of pr(w0) */
    int    *resultcode /* 0 for successful completion */
)

{
    int n,k,m,c,i,x;
    int gi;
    int wxi;
    gfnptr gfn;
    double asum = 0;
    double *gk0;
    double lambda;
    double par[4];
    double stdint;
    int cumk[1001];
    int nk = 0;
    double *pmix;      /* proportion in each mixture */
    int gpar = 2;

    *resultcode = 1;                   /* generic failure */

    /* determine number of polygons or transects */
    /* for polygon or transects detectors, kk is vector ending in zero */
    if ((*detect==6) | (*detect==7)) {
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
    }
    else nk = *kk;

    /* Allocate space for array of naive detection probability */
    gk0 = (double *) R_alloc(*cc0 * nk * *mm, sizeof (double));
    pmix = (double *)  R_alloc(*nc * *nmix, sizeof (double));
    for (i=0; i< *nc * *nmix; i++) pmix[i] = 1; /* default */

    /* 
        *fn may take values -
        0  halfnormal
        1  hazard rate
        2  exponential
        3  compound halfnormal
        5  w-exponential
        10 signal strength (signal detectors only)
        11 binary signal strength
    */

    if (*fn == 0)       { gfn = ghn; }
    else if (*fn == 1)  { gfn = ghz; gpar++; }
    else if (*fn == 2)  { gfn = ghe; }
    else if (*fn == 3)  { gfn = ghnc; gpar++; }
    else if (*fn == 5)  { gfn = ghf; gpar++; }
    else if (*fn == 9)  { gfn = gsigbin; }
    else if (*fn == 10) { gfn = gsig; gpar++; }
    else if (*fn == 11) { gfn = gsigsph; gpar++; }
    else return;      /* return if detection function not recognised */
    if (*nmix>1) gpar++;

    if (*detect == 0) {
        for (c=0; c<*cc0; c++) {
            for (k=0; k<*kk; k++) {
                for (m=0; m<*mm; m++) {        
                        gi = i3(c, k, m, *cc0, nk);
                        gk0[gi] = gfn(k, m, c, gsb0val, *cc0, traps, mask, nk, *mm, *cut);
                }
            }
        }
    }
    else if (((*detect >= 1) && (*detect <= 5)) | (*detect==8)) {
        for (c=0; c<*cc0; c++) {
            for (k=0; k<*kk; k++) {
                for (m=0; m<*mm; m++) {
                    lambda = gfn(k, m, c, gsb0val, *cc0, traps, mask, nk, *mm, *cut);
                    gi = i3(c, k, m, *cc0, nk);
                    if (*binomN == 0)
                        gk0[gi] = 1 - gpois (0, lambda, 0);
                    else if (*binomN == 1)
                        gk0[gi] = lambda;
                    else if (*binomN < 0)
                        gk0[gi] = 1 - gnbinom (0, *binomN, lambda, 0);
                    else 
                        gk0[gi] = 1 - gbinom (0, *binomN, lambda / *binomN, 0);
                }
            }
        }
    }
    else if (*detect == 6) {
        for (c=0; c<*cc0; c++) {
            par[0] = 1;
            par[1] = gsb0val[*cc0 + c];
            par[2] = gsb0val[2* *cc0 + c];   
            stdint = gintegral(*fn, par);
            for (k=0; k<nk; k++) {                 /* over polygons */
                for (m=0; m<*mm; m++) {
                    lambda = integral2D (*fn, m, c, gsb0val, *cc0, traps, 
                        mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
                    gi = i3(c, k, m, *cc0, nk);
                    gk0[gi] = 1 - gpois (0, lambda, 0);
                }
            }
        }
    }
    else if (*detect == 7) {
        for (c=0; c<*cc0; c++) {
            par[0] = gsb0val[c];
            par[1] = gsb0val[*cc0 + c];
            par[2] = gsb0val[2* *cc0 + c];   
            for (k=0; k<nk; k++) {                 /* over transects */
                for (m=0; m<*mm; m++) {
                    lambda = integral1D (*fn, m, c, gsb0val, *cc0, traps, 
                        mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm);
                    gi = i3(c, k, m, *cc0, nk);
                    gk0[gi] = 1 - gpois (0, lambda, 0);
                }
            }
        }
    }
    else error ("unrecognised detector type in external C fn integralprw1");

    if (*nmix>1) {
        for (n=0; n<*nc; n++) {          
            for (x=0; x<*nmix; x++) {        
                wxi = i4(n,0,0,x,*ncol,*ss,*kk);
                c = gsb0[wxi] - 1; 
                pmix[*nmix * n + x] = gsb0val[*cc0 * (gpar-1) + c];
           }        
       }
    }

    for (n=0; n<*nc; n++) {            /* CH numbered 0 <= n < *nc */
        if ((*ncol > 1) | (n == 0)) {  /* no need to repeat if constant */     
            if ((n+1) > *ncol) {       /* groups not implemented */
                *resultcode = 3;  
                return; 
            }
            asum = 0;
            for (x=0; x<*nmix; x++) {
                for (m=0; m<*mm; m++) {
                    asum += pmix[*nmix * n + x] * pndot (m, n, 1, *ss, x, *ncol, gsb0, gk0, *ss, nk, *cc0, *nmix);   
                }
            }
        }
        a[n] = *area * asum;
    }
    *resultcode = 0;                   /* successful completion */
}

/*==============================================================================*/

double hxy (
    int c, 
    int m,
    double gk[],
    int cc,
    int kk)
/*
    Sum hazard over traps for animal at m
    for parameter combination c
*/
{
    int k;
    double pp;
    double sumhz = 0;
    for (k=0; k<kk; k++) {
        /* probability caught in trap k 0 <= k < kk */
        pp = gk[i3(c, k, m, cc, kk)];
        if (pp > (1-fuzz))             /* pp close to 1.0 - approx limit */
            pp = huge;                 /* g0 very large (effectively infinite hazard) */
        else {
            if (pp <= 0) pp = 0;
            else pp = -log(1-pp);
        }
        sumhz += pp;                   /* accumulate hazard over traps */
    }
    return (sumhz);
}
/*==============================================================================*/

double hxys (
    int c, 
    int m, 
    int s, 
    int gsb[], 
    double gk[],
    int cc,
    int nc, 
    int ss, 
    int kk
)
/*
    Sum hazard over traps for animal at m
    for parameter combination c
    with traps in use on occasion s
*/
{
    int k;
    double pp;
    double sumhz = 0;
    for (k=0; k<kk; k++) {
        /* check animal 0 to see if this trap set on occ s */
        /* skip if not set */
        if (gsb[i4(0, s, k, 0, nc, ss, kk)] > 0) {
            /* probability caught in trap k 0 <= k < kk */
            pp = gk[i3(c, k, m, cc, kk)];
            if (pp > (1-fuzz))         /* pp close to 1.0 - approx limit */
                pp = huge;             /* g0 very large (effectively infinite hazard) */
            else {
                if (pp <= 0) pp = 0;
                else pp = -log(1-pp);
            }
            sumhz += pp;               /* accumulate hazard over traps */
        }
    }
    return (sumhz);
}
/*========================================================*/

double prwimulti 
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[], 
    double gk[], int binomN, double hxytemps[], int cc, int nc, int kk, int ss, int mm, int nmix, 
    gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < *nc)
    given that animal's range centre is at m
    MULTI-CATCH DETECTOR
*/
{
    int s;                             /* index of occasion  0 <= s < *ss */
    int k;                             /* index of trap      0 <= k < *kk */
    int c;
    int gi;
    int dead = 0;
    double htemp;
    double pks;
    double result = 1.0;
  
    for (s=s1-1; s<s2; s++)
    {
        k = w[nc * s + n];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        if (k > 0) {
            c = gsb[i4(n, s, k-1, x, nc, ss, kk)] - 1;  
            htemp = hxytemps[i3(c, m, s, cc, mm)];
            if (htemp < fuzz) { result = 0; break; }
            gi = i3(c, k-1, m, cc, kk);
            pks = -log (1 - gk[gi]);
            pks *= (1-expmin(-htemp)) / htemp;   
        }
        else  {
            c = gsb[i4(n, s, 0, x, nc, ss, kk)] - 1;   /* use combination for k=0 - as good as any */
            htemp = hxytemps[i3(c, m, s, cc, mm)]; 
            if (htemp < fuzz) { result = 0; break; }
            pks = expmin(-htemp);      /* Not captured */
        }
        result *= pks;
        if (dead) break;
    }
    return (result);
}  
/*=============================================================*/

double prwicount 
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[], 
     double gk[], int binomN, double detspec[], int cc, int nc, int kk, int ss, int mm, int nmix,
     gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    COUNT DETECTOR, includes PROXIMITY
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of trap      0 <= k < kk  */
    int c, gi, wi, wxi;
    int count;
    int dead = 0;
    double result = 1.0;  
    
    for (s=s1-1; s<s2; s++) { 
        for (k=0; k<kk; k++) {
            wi = i3(n, s, k, nc, ss);
            wxi = i4(n, s, k, x, nc, ss, kk);
            count = w[wi];    
            if (count<0) {dead=1; count=-count;}
            c = gsb[wxi] - 1;  
            if (c >= 0) {                                       /* skip if this trap not set */
                gi  = i3(c, k, m, cc, kk);
                if (binomN == 0)                                /* Poisson */
                    result *= gpois(count, gk[gi], 0);        
                else if (binomN == 1) {                         /* Bernoulli */
                    if (count>0)
                        result *= gk[gi];                     
                    else
                        result *= 1 - gk[gi];        
                }
                else if (binomN < 0)                           /* negative binomial */
                    result *= gnbinom (count, binomN, gk[gi], 0);
                else                                            /* binomial */ 
                    result *= gbinom(count, binomN, gk[gi]/binomN, 0);        
                if (result < minp) {result = minp; break;}      /* truncate */
            }
        }
        if (result <= minp) {result = minp; break;}             /* truncate */
        if (dead==1) break;
    }
    return (result);
}
/*=============================================================*/

double prwisignal     
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[], 
     double gk[], int binomN, double detspec[], int cc, int nc, int kk, int ss, int mm, int nmix,
     gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < nc) 
    given that animal's range centre is at m    
    SIGNAL STRENGTH DETECTOR                    
*/

{
    int s;   /* index of occasion  0 <= s < ss */
    int k;   /* index of trap      0 <= k < kk */
    int j;   /* index of detection */ 
    int c, wi, wxi, gi;
    double result = 1.0;      
    double mu, gam, sdS, y;
    int start = 0;
    int count = 0;
    double cut;
    int spherical;

    cut = detspec[0];
    spherical = round (detspec[1]);
    for (s=s1-1; s<s2; s++) { 
        for (k=0; k<kk; k++) {           
            wxi = i4(n, s, k, x, nc, ss, kk);
            c = gsb[wxi] - 1;  
            if (c >= 0) {                                           /* skip if this trap not set */
                wi = i3(n, s, k, nc, ss);
                count = abs(w[wi]);                                 /* ignore 'deaths */
                if (count==0) {
                    gi = i3(c,k,m,cc,kk);
                    if (binomN == 0)                                /* Poisson */
                        result *= gpois(0, gk[gi], 0);        
                    else if (binomN == 1)                           /* Bernoulli */
                        result *= 1 - gk[gi];
                    else if (binomN < 0)                            /* negative binomial */
                        result *= gnbinom (0, binomN, gk[gi], 0);
                    else                                            /* binomial */ 
                        result *= gbinom(0, binomN, gk[gi]/binomN, 0);        
                }
                else {
                    start = round(detspec[wi+2]);               
                    for (j=0; j<count; j++) {
                        if (spherical == 0) {
                            mu = mufn (k, m, gsbval[c], gsbval[cc + c],
                                traps, mask, kk, mm); 
                        }
                        else {
                            mu = mufnsph (k, m, gsbval[c], gsbval[cc + c],
                                traps, mask, kk, mm); 
                        }
                        sdS = gsbval[cc * 2 + c];
                        gam = (cut - mu) / sdS; 
                        y   = signal[start+j] - cut + gam * sdS;
                        result *= dnorm(y, 0, sdS, 0);                  /* signal strength */
                    }
                }
                if (result < minp) {result = minp; break;}          /* truncate */
            }
        }
        if (result <= minp) {result = minp; break;}                 /* truncate */
    }
    return (result);
}
/*=============================================================*/

double prwitimes     
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[], 
     double gk[], int binomN, double detspec[], int cc, int nc, int kk, int ss, int mm, int nmix, 
     gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < nc) 
    given that animal's range centre is at m    
    TIMES DETECTOR  (uses 'signal' vector for times)
*/

/* R-exts.pdf "the exponential and gamma distributions are parametrized 
   by scale rather than rate" see code in src\nmath\... */

{
    int s;   /* index of occasion  0 <= s < ss */
    int k;   /* index of trap      0 <= k < kk */
    int c, wi, wxi, gi;
    double lambda;
    double result = 1.0;      
    int count = 0;

    /*  for times ... */
    int j; 
    int start = 0;
    double time0;
    double y = 0;

    for (s=s1-1; s<s2; s++) { 
        for (k=0; k<kk; k++) {           
            wxi = i4(n,s,k,x,nc,ss,kk);
            wi = i3(n,s,k,nc,ss);
            c = gsb[wxi] - 1;  
            if (c >= 0) {                                       /* skip if this trap not set */
                count = abs(w[wi]);
                gi = i3(c,k,m,cc,kk);
                lambda = gk[gi];

/*
    The following 'intuitive' code attempts to implement a continuous-time
    model. This appears to be ill-conceived because in a homogeneous
    model the times contain no additional information (e.g. Nayak 1988 
    Biometrika 75:113). 

    It is unclear whether the final likelihood component is needed.
    A small trial without it suggests the asymptotic variances estimates
    are slightly graeter than the corresponding 'count' detector estimates.

    'times' detectors are therefore not mentioned in current documentation.

    MGE 2009 12 07
*/
                if (binomN == 0)                              
                    result *= gpois(count, lambda, 0);        
                else if (binomN == 1) {
                    if (count==0)      
                        result *= 1 - lambda;
                    else
                        result *= lambda;
                }
                else if (binomN < 0)                           
                    result *= gnbinom (count, binomN, lambda, 0);
                else                                            
                    result *= gbinom(count, binomN, lambda/binomN, 0);        

                if (count>0) {
                    start = round(detspec[wi]);
                    time0 = 0.0;
                    for (j=0; j<count; j++) {
                        y = signal[start+j]-time0;
                        result *= dexp(y, 1/lambda, 0) / pexp (1-time0, 1/lambda, 1, 0);
                        time0 = signal[start+j];
                    }
                    /* result *= 1- pexp (1-time0, 1/lambda, 1, 0);   no detection after last */
                } 
                if (result < minp) {result = minp; break;}      /* truncate */
            }
        }
        if (result <= minp) {result = minp; break;}             /* truncate */
    }
    return (result);
}
/*=============================================================*/

double prwipolygon 
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[], 
     double gk[], int binomN, double detspec[], int cc, int nc, int kk, int ss, int mm, int nmix,
     gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)

/*
    Likelihood component due to capture history n (0 <= n < nc)
    given that animal's range centre is at m
    POLYGON OR TRANSECT DETECTOR  
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of part 0 <= k < nk  */
    int j;   /* index of xy record */
    int c, wxi, wi, gi;
    long count;
    int dead = 0;
    double result = 1.0; 
    double geval; 
    int nk;
    int nd;
    int start;

    nk = round(detspec[0]);
    nd = round(detspec[1]);
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<nk; k++) {
            wi = i3(n,s,k,nc,ss);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            wxi = i4(n,s,k,x,nc,ss,kk);
            c = gsb[wxi] - 1;  
            if (c >= 0) {                                       /* skip if this polygon not used */
                gi  = i3(c,k,m,cc,nk);     
                result *= gpois(count, gk[gi], 0);              /* Poisson - lambda is cues/animal x cue detection */
                if (result > minp) {                            /* avoid underflow 2009 11 16 */             
                    start = round(detspec[wi+cc+2]);
                    for (j=start; j < start+count; j++) {
                        geval = gfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0);
                        result *= geval / (gk[gi] * detspec[2+c]);          
                    }
                }
            }
        }
        if (result <= minp) {result = minp; break;}  
        if (dead==1) break;
    }
    return (result);
}
/*=============================================================*/

void pdotpoly (double *xy, int *nxy, double *traps, int *nk, 
    int *kk, int *fn, int *nocc, double *par, double *value)
{
    double stdintegral;
    int k;
    int i;
    double lambda;
    double sumlambda;
    int *cumk;
    cumk = (int *) R_alloc(*nk+1, sizeof(int));
    cumk[0] = 0;
    for (i=0; i<*nk; i++) 
        cumk[i+1] = cumk[i] + kk[i];
    stdintegral = gintegral(*fn, par);
    for (i=0; i<*nxy; i++) {
        sumlambda = 0.0;
        for (k=0; k<*nk; k++) {                        /* over parts */
            lambda = integral2D (*fn, i, 0, par, 1, traps, xy, 
                cumk[k], cumk[k+1]-1, cumk[*nk], *nxy) / stdintegral;
            sumlambda += lambda;            
        }   
        value[i] = 1 - pow(gpois (0, sumlambda, 0), *nocc);
    }
}
/*=============================================================*/

void pdottransect (double *xy, int *nxy, double *traps, int *nk, 
    int *kk, int *fn, int *nocc, double *par, double *value)
{
    int k;
    int i;
    double lambda;
    double sumlambda;
    int *cumk;
    cumk = (int *) R_alloc(*nk+1, sizeof(int));
    cumk[0] = 0;
    for (i=0; i<*nk; i++) 
        cumk[i+1] = cumk[i] + kk[i];
    for (i=0; i<*nxy; i++) {
        sumlambda = 0.0;
        for (k=0; k<*nk; k++) {                        /* over transects */
            lambda = integral1D (*fn, i, 0, par, 1, traps, xy, 
                cumk[k], cumk[k+1]-1, cumk[*nk], *nxy);
            sumlambda += lambda;            
        }   
        value[i] = 1 - pow(gpois (0, sumlambda, 0), *nocc);
    }
}
/*=============================================================*/

void secrloglik (
    int    *like,        /* likelihood 0 full, 1 conditional */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
    int    *distrib,     /* distribution of nc 0 Poisson, 1 binomial */
    int    *w,           /* capture histories (1:nc, 1:s, 1:k) */
    double *xy,          /* xy coordinates of polygon records */
    double *signal,      /* signal strength vector, or times */
    int    *grp,         /* group number for 0<=n<*nc   [full likelihood only] */
    int    *nc,          /* number of capture histories */
    int    *ss,          /* number of occasions */
    int    *kk,          /* number of traps */
    int    *mm,          /* number of points on mask */
    int    *gg,          /* number of groups */
    int    *nmix,        /* number of mixtures */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *Dmask,       /* density at each point on mask, possibly x group */
    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    double *gsb0val,     /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *cc0,         /* number of g0/sigma/b combinations for naive animals */
    int    *gsb,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *gsb0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive animal] */
    double *area,        /* area associated with each mask point (ha) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *cut,         /* transformed signal strength threshold for detection */
    double *minprob,     /* minimum value of P(detection history) */
    double *a,           /* a(theta) */
    double *value,       /* return value integral of pr(w0) */
    int    *resultcode   /* 0 if OK */
)
{
    int    i,n,g,k,m,c,s,x;
    int    ng[*gg];       /* number per group */
    int    nSK;
    int    wi;
    int    wxi;
    int    gi;
    int    notset = 0;
    int    indiv = 0;     /* indicator for whether detection varies between individuals (CL only) */
    double asum[maxnmix];      
    double *pmix;
    double temp, tempsum;
    double sumD[*gg];     /* 2009 08 21  for binomial */
    double sumDp[*gg];

    double *gk;
    double *gk0;
    double *detspec;
    int    *start = NULL;
    int    nval;
    double lambda;

    gfnptr gfn;
    prwfnptr prwfn;
    int cumk[1001];
    double par[4];
    double stdint;
    double templog;
    int nk = 0;
    int nd = 0;
    int gpar = 2;    /* number of 'g' (detection) parameters */

    /*===============================================================*/

    /* MAINLINE */

    *resultcode = 1;  /* generic failure code */
                      /* reset to 0 at end */

    /* determine number of polygons if polygon detector */
    /* for polygon detectors, kk is vector ending in zero */
    if ((*detect==6) | (*detect==7)) {
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
    }
    else 
        nk = *kk;

    if ((*detect>=5) && (*detect<=8)) {
        /* start[z] indexes the first row in xy (or element in signal) 
           for each possible count z, where z is w-order (isk) */
        start = (int *) R_alloc(*nc * *ss * nk, sizeof(int));
        for (k=0; k<nk; k++) {
            for (s=0; s< *ss; s++) {
                for (i=0; i< *nc; i++) {
                    wi = i3(i,s,k,*nc,*ss);
                    start[wi] = nd;
                    nd += abs(w[wi]);
                } 
            }
        }
    }

    if (*detect==0) 
        nval = *cc * *mm * *ss;
    else if (*detect==5) 
        nval = 2 + *nc * *ss * nk * *nmix;
    else if ((*detect==6) | (*detect==7))
        nval = 2 + *cc + *nc * *ss * nk * *nmix;
    else if (*detect==8) 
        nval = *nc * *ss * nk * *nmix;
    else 
        nval = 4;    /* 1-4, mostly not used */

    /* use R_alloc for robust exit on interrupt */
    gk = (double *) R_alloc(*cc * nk * *mm, sizeof(double));
    gk0 = (double *) R_alloc(*cc0 * nk * *mm, sizeof(double));
    detspec = (double *) R_alloc(nval, sizeof(double));
    pmix = (double *)  R_alloc(*nc * *nmix, sizeof (double));

    /* 
        *fn may take values -
        0  halfnormal
        1  hazard rate
        2  exponential
        3  compound halfnormal
        5  w-exponential
        10 signal strength (signal detectors only)
        11 binary signal strength
    */

    if (*fn == 0)       { gfn = ghn; }
    else if (*fn == 1)  { gfn = ghz; gpar++; }
    else if (*fn == 2)  { gfn = ghe; }
    else if (*fn == 3)  { gfn = ghnc; gpar++; }
    else if (*fn == 5)  { gfn = ghf; gpar++; }
    else if (*fn == 9)  { gfn = gsigbin; }
    else if (*fn == 10) { gfn = gsig; gpar++; }
    else if (*fn == 11) { gfn = gsigsph; gpar++; }
    else return;
    if (*nmix>1) gpar++;

    /*===============================================================*/
    /* populate pre-computed arrays */

    /*
        *detect may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  binary area search
        4  count  area search
        5  signal detectors
        6  polygon detector
        7  transect detector
        8  times 
    */

    /* dedicated detectors ... */
    if (*detect==1) *binomN = 1;  /* Bernoulli */
    if (*detect==3) *binomN = 1;  /* Bernoulli */

    if (*detect == 0) {          
        for (k=0; k<*kk; k++) {
            for (m=0; m<*mm; m++) {
                for (c=0; c<*cc0; c++) {
                    gi = i3(c,k,m,*cc0,nk);
                    gk0[gi] = gfn(k, m, c, gsb0val,
                        *cc0, traps, mask, nk, *mm, *cut);
                }
                for (c=0; c<*cc; c++) {
                    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = gfn(k, m, c, gsbval,
                        *cc, traps, mask, nk, *mm, *cut);
                }
            }
        }
    }
    else if (((*detect >= 1) && (*detect <= 5)) | (*detect==8)) {
        for (k=0; k<*kk; k++) {
            for (m=0; m<*mm; m++) {
                for (c=0; c<*cc0; c++) {
                    lambda = gfn(k, m, c, gsb0val, *cc0, traps, mask, *kk, 
                        *mm, *cut);
                    gi = i3(c,k,m,*cc0,nk);
                    if (*binomN == 0)
                        gk0[gi] = 1 - gpois (0, lambda, 0);
                    else if (*binomN == 1)
                        gk0[gi] = lambda;
                    else if (*binomN < 0)
                        gk0[gi] = 1 - gnbinom (0, *binomN, 
                            lambda, 0);
                    else  /* (*binomN > 1) */
                        gk0[gi] = 1 - gbinom (0, *binomN, lambda / *binomN, 0);
                }
                for (c=0; c<*cc; c++) {
                    gi = i3(c,k,m,*cc,nk); 
                    gk[gi] = gfn(k, m, c, gsbval, *cc, 
                        traps, mask, nk, *mm, *cut);
                } 
            }
        }
    }
    else if (*detect == 6) {
        for (c=0; c<*cc0; c++) {
            par[0] = 1;
            par[1] = gsb0val[*cc0 + c];
            par[2] = gsb0val[2* *cc0 + c];   
            stdint = gintegral(*fn, par);
            for (k=0; k<nk; k++) {               /* over parts */
                for (m=0; m<*mm; m++) {
                    lambda = integral2D (*fn, m, c, gsb0val, *cc0, traps, mask, 
                        cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint ;
                    gi = i3(c,k,m,*cc0,nk);
                    gk0[gi] = 1 - gpois (0, lambda, 0);
                }
            }
        }
        R_CheckUserInterrupt();
        for (c=0; c<*cc; c++) {
            par[0] = 1;
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];   
            detspec[2+c] = gintegral(*fn, par);  /* passed to prwipolygon */
            for (k=0; k<nk; k++) {               /* over parts */
                for (m=0; m<*mm; m++) {
                    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = integral2D (*fn, m, c, gsbval, *cc, 
                        traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / detspec[2+c];                    
                }
            }
        }
    }
    else if (*detect == 7) {
        for (c=0; c<*cc0; c++) {
            for (k=0; k<nk; k++) {               /* over transects */
                for (m=0; m<*mm; m++) {
                    lambda = integral1D (*fn, m, c, gsb0val, *cc0, traps, mask, 
                        cumk[k], cumk[k+1]-1, cumk[nk], *mm);
                    gi = i3(c,k,m,*cc0,nk);
                    gk0[gi] = 1 - gpois (0, lambda, 0);
                }
            }
        }
        R_CheckUserInterrupt();
        for (c=0; c<*cc; c++) {
            detspec[2+c] = 1;                    /* passed to prwipolygon */
            for (k=0; k<nk; k++) {               /* over transects */
                for (m=0; m<*mm; m++) {
                    gi = i3(c,k,m,*cc0,nk);
                    gk[gi] = integral1D (*fn, m, c, gsbval, *cc, 
                        traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm);                    
                }
            }
        }
    }
    else error ("unrecognised detector type in external C fn secrloglik");

    R_CheckUserInterrupt();

    /*===============================================================*/

    prwfn = prwicount;   /* default */

    if (*detect == 0) {             
        prwfn = prwimulti; 
        /* check if any traps not set */
        /* and decide whether to use full or collapsed hxy */
        nSK = *nc * *ss * *kk;
        for (i=0; i<nSK; i++) 
            if (gsb[i]==0) notset++;
        if (notset==0) {
            for (c=0; c<*cc; c++) {
                for (m=0; m<*mm; m++) { 
                    detspec[i3(c,m,0,*cc,*mm)] = hxy(c,m,gk, *cc, *kk);
                    for (s=1; s<*ss; s++) {
                        detspec[i3(c,m,s,*cc,*mm)] = detspec[i3(c,m,0,*cc,*mm)];
                    } 
                }
            }
        }
        else {    /* some traps not set: occasion-specific total hazard */
          for (s=0; s<*ss; s++)  
            for (c=0; c<*cc; c++)
              for (m=0; m<*mm; m++) 
                  detspec[i3(c,m,s,*cc,*mm)] = hxys(c,m,s, gsb, gk, *cc, *nc, *ss, *kk);
        }
    }
    else if (*detect == 1) { prwfn = prwicount; }    /* Bernoulli */
    else if (*detect == 2) { prwfn = prwicount; }    /* Poisson (0) or binomial or neg binomial */     
    else if (*detect == 3) { prwfn = prwicount; }    /* Bernoulli */
    else if (*detect == 4) { prwfn = prwicount; }    /* Poisson (0) or binomial or neg binomial */     
    else if (*detect == 5) { 
        prwfn = prwisignal; 
        detspec[0]= *cut; 
        detspec[1]= (*fn == 11);     /* spherical */ 
        for (i=0; i< (*nc* *ss * nk); i++) 
            detspec[2+i] = (double) start[i];
    }
    else if ((*detect == 6) | (*detect == 7)) { 
        prwfn = prwipolygon;
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        for (i=0; i< (*nc* *ss * nk); i++) 
            detspec[2+*cc+i] = (double) start[i];
    }
    else if (*detect == 8) { 
        prwfn = prwitimes; 
        for (i=0; i< (*nc* *ss * nk); i++) 
            detspec[i] = (double) start[i];
    }
    for (i=0; i < *nc * *nmix; i++) pmix[i] = 1; /* default */
    if (*nmix>1) {
        for (n=0; n<*nc; n++) {          
            for (x=0; x<*nmix; x++) {        
                wxi = i4(n,0,0,x,*nc,*ss,nk);
                c = gsb[wxi] - 1; 
                if (*like==0) {
                    g = grp[n]-1;
                    pmix[*nmix * g + x] = gsbval[*cc * (gpar-1) + c];
                }
                else
                    pmix[*nmix * n + x] = gsbval[*cc * (gpar-1) + c];
            }        
        }
    }
    /*===============================================================*/
    if (*like==1)    /* Conditional likelihood */
    {
        /*
           check if we need to consider variation among individuals
           i.e. check if detection parameters constant for given s,k
        */
        indiv = 0;
        for (s=0; s<*ss; s++) { 
            for (k=0; k<nk; k++) {   
                for (x=0; x<*nmix; x++) {
                    wxi = i4(0,s,k,x,*nc,*ss,*kk);
                    i = gsb0[wxi];
                    for (n=1; n<*nc; n++) {
                        wxi = i4(n,s,k,x,*nc,*ss,*kk);
                        if (i != gsb0[wxi]) {
                            indiv = 1; break;
                        }
                    }
                }
            }
        }

        if (indiv == 0) {
            /* save time by doing this once, rather than inside n loop */
            for (x=0; x<*nmix; x++) {        
                asum[x] = 0;
                for (m=0; m<*mm; m++) {
                    asum[x] += pndot (m, 0, 1, *ss, x, *nc, gsb0, gk0, *ss, nk, *cc0, *nmix);   /* all individuals the same */
                }
            }
        }
        /* else asum calculated for each individual in loop below */

        *value = 0;
        /* Loop over individuals... */
        for (n=0; n<*nc; n++) {                      /* CH numbered 0 <= n < *nc */
            a[n] = 0; 
            tempsum = 0;
            for (x=0; x<*nmix; x++) {
                temp = 0;
                if (indiv > 0) 
                    asum[x] = 0;
                for (m=0; m<*mm; m++) {
                    if (indiv > 0) 
                        asum[x] += pndot (m, n, 1, *ss, x, *nc, gsb0, gk0, *ss, nk, *cc0, *nmix);  
                    temp += prwfn (m, n, 1, *ss, x, w, xy, signal, gsb, gk, *binomN, detspec, 
                        *cc, *nc, nk, *ss, *mm, *nmix, gfn, gsbval, traps, mask, *minprob);
                }
                a[n] += pmix[*nmix * n + x] * asum[x];    
                tempsum += pmix[*nmix * n + x] * temp;
            }    /* end of loop over mixtures */

            templog = log(tempsum/a[n]);
            a[n] = *area * a[n];
           
            if (R_FINITE(templog)) 
                *value += templog;
            else {
                *resultcode = 9;
                return; 
            }              
            R_CheckUserInterrupt();
        }        /* end of loop over individuals */
    }
    /*-------------------------------------------------------------------------------------------*/

    else {  /* *like==0  Full likelihood */

        for (g=0; g<*gg; g++) {
            ng[g] = 0;  /* zero for later */
            sumD[g] = 0;
            sumDp[g] = 0;
            for (m=0; m<*mm; m++)  {
                sumD[g] += Dmask[*mm * g + m];
                for (x=0; x<*nmix; x++)
                    sumDp[g] += pmix[*nmix * g + x] * pndot (m, g, 1, *ss, x, *gg, gsb0, gk0, 
                        *ss, nk, *cc0, *nmix) * Dmask[*mm * g + m];
            }
        }
        *value = 0;
        for (n=0; n<*nc; n++) {  /* CH are numbered 0 <= n < *nc in C code */
            g = grp[n]-1;
            ng[g]++; 
            temp = 0;
            for (x=0; x<*nmix; x++) {
                for (m=0; m<*mm; m++) {
                    temp += pmix[*nmix * g + x] * prwfn (m, n, 1, *ss, x, w, xy, signal, gsb, 
                        gk, *binomN, detspec, *cc, *nc, nk, *ss, *mm, *nmix, gfn, 
                        gsbval, traps, mask, *minprob) * Dmask[*mm * g + m];    
                }
            }
            templog = log(temp);
            if (R_FINITE(templog)) 
                *value += templog;
            else {
                *resultcode = 9;
                return; 
            }
            R_CheckUserInterrupt();
        }

        for (g=0; g<*gg; g++) {
            *value -= ng[g] * log(sumDp[g]);
            /* Poisson */
            if (*distrib==0) *value += gpois(ng[g], sumDp[g] * *area, 1);
            /* binomial */
            /* note must round 'size' argument to an integer */
            /* note restriction to integer values causes asymptotic variance of density to fail */
            if (*distrib==1) *value += gbinom (ng[g], round(sumD[g] * *area), sumDp[g]
                / sumD[g], 1);                                       
        }
    }

    *resultcode = 0;   /* successful termination secrloglik */
}
/*==============================================================================*/

/*
    trappingXXXX routines perform simulated sampling of 2D popn with various 
    detector types
*/

double pfn (
    int fn,
    double d2val, 
    double g0,
    double sigma,
    double z,
    double area)
{
    double p = -1;
    double gam = 0;
    if (fn == 0) p = g0 * expmin(-d2val / 2 / sigma / sigma);
    else if (fn == 1) p = g0 * (1 - expmin(- pow (sqrt(d2val) / sigma, - z)));
    else if (fn == 2) p = g0 * expmin(- sqrt(d2val) / sigma);
    else if (fn == 3) { if (sqrt(d2val) <= sigma) p = g0; else p = 0; }
    else if (fn == 5) { if (sqrt(d2val) <= z) p = g0; else p = g0 * expmin(- (sqrt(d2val)-z) / sigma); }
    else if (fn == 11) { gam = -(g0 + sigma * sqrt(d2val)); p = pnorm(gam,0,1,0,0); }
    else error ("requested detection function not valid in this context");
    return (p);
}

double randomtime (double p)
/* return random event time for event with probability p */
{
    double minprob = 1e-5;
    double lambda;
    double random_U;

    if (p < minprob)
        return(huge);                        /* ignore trivial p/lambda */
    else if (p >= 1.0)
        return (-Random());                  /* trick to spread P=1 */
    else {
        lambda   = -log(1-p);                /* rate parameter */
        random_U = Random();
        if (random_U <= 0)                   /* trap for zero */
            return(huge);
        else 
            return (-log(random_U)/lambda);   /* random exponential e.g. Ripley 1987 Algorithm 3.2 p 55 */
    }
}
/*==============================================================================*/

    void probsort (int n, struct trap_animal tran[])
    /*
        Sort using Shell algorithm see Press et al 1989 p 257 
        tran is an array of trap_animal records 
    */
    {
       double aln2i = 1.442695022;
       double tiny  = 1.0e-5;
       int nn,m,lognb2,l,k,j,i;
       struct trap_animal t; 
    
       lognb2 = trunc(log(n)*aln2i+tiny);
       m = n;
       for (nn=1; nn<=lognb2; nn++)
       {
          m = m / 2;
          k = n-m;
          for (j=1; j<=k; j++)
          { 
             i = j;
    lab1:    l = i+m;
             if (tran[l-1].time < tran[i-1].time) 
             {
                t = tran[i-1];
                tran[i-1] = tran[l-1];
                tran[l-1] = t;
                i = i-m;
                if (i >= 1)  goto lab1;
             }
          }
       }
    }    /* end of probsort */
/*==============================================================================*/

void trappingcount (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    int    *used,      /* ss x kk array of 0/1 codes for usage */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,    /* 0 poisson, 1 Bernoulli, or number of trials for 'count' 
                          detector modelled with binomial */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)

{

    double d2val;
    double theta;
    int    i,j,k,l,s;
    int    nc;
    int    count;

    *resultcode = 1;
    nc = 0;
    GetRNGstate();

    for (i=0; i<*N; i++) caught[i] = 0; 
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            for (k=0; k<*kk; k++) {
                if (used[s * *kk + k]) {                        /* 2009 11 09 */            
                    d2val = d2(i,k, animals, traps, *N, *kk);
                    theta = pfn(*fn, d2val, g0[s], sigma[s], z[s], 0);
                    if (theta>0) {
                        if (*binomN == 0)
                            count = rpois(theta);
                        else if (*binomN == 1) {
                            if (Random() < theta)
                                count = 1;
                            else
                                count = 0; 
                        }
                        else if (*binomN < 0) {
                            /* must use 'size, prob' parameters */
                            /* prob = size / (size + mu) */
                            *binomN = abs(*binomN);
                            count = rnbinom(*binomN, *binomN / (*binomN+theta));   
                        }
                        else
                            count = rbinom(*binomN, theta / *binomN);   
                        if (count>0)
                        {
                             if (caught[i]==0)                  /* first capture of this animal */
                             {
                                 nc++; 
                                 caught[i] = nc;
                                 for (j=0; j<*ss; j++)
                                   for (l=0; l<*kk; l++)
                                     value[*ss * ((nc-1) * *kk + l) + j] = 0;
                             } 
                             value[*ss * ((caught[i]-1) * *kk + k) + s] = count;
                        }
                    }
                }
            }
        }     
    }
    *n = nc;

    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingpolygon (
    double *lambda,      /* Parameter : NOT expected detection events per hectare */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *npoly,       /* number of different polygons */
    int    *kk,          /* number of vertices + 1 (assumed closed) for each polygon */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *maxone,      /* maximum of one detection per animal per occasion */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)

{
    int    i,j,k,l,s,t;
    int    nc = 0;
    int    np = 1;       /* number of points each call of gxy */
    int    nd = 0;
    int    count;
    double par[3];
    double w;
    int    maxdet;
    int    g=0;
    int    *gotcha;
    double xy[2]; 
    gotcha = &g;
    int cumk[maxnpoly+1];   /* limit maxnpoly polygons */
    int sumk;
    int n1,n2;
    double *workXY;
    int    *sortorder;
    double *sortkey;

    *resultcode = 1; 
    GetRNGstate();
    if (*npoly>maxnpoly) {
        *resultcode = 2;
        return;
    }
    cumk[0] = 0;
    for (k =0; k<*npoly; k++) cumk[k+1] = cumk[k] + kk[k];
    sumk = cumk[*npoly];
    maxdet = *N * *ss * *npoly * maxperpoly;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));

    for (i=0; i<*N; i++) caught[i] = 0;                           
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            if (*maxone)
                count = (Random() < lambda[s]);
            else
                count = rpois(lambda[s]);              /* number of cues */
            w = 5 * sigma[s];
            par[0] = 1;
            par[1] = sigma[s];
            par[2] = z[s];
            for (j=0; j<count; j++) {             
                gxy (&np, fn, par, &w, xy);            /* simulate location */
                xy[0] = xy[0] + animals[i];
                xy[1] = xy[1] + animals[*N + i];
                for (k=0; k<*npoly; k++) {             /* each polygon */
                    n1 = cumk[k]; 
                    n2 = cumk[k+1]-1;
                    inside(xy, &n1, &n2, &sumk, traps, gotcha);  /* assume closed */
                    if (*gotcha > 0) {
                        if (caught[i]==0) {            /* first capture of this animal */
                            nc++; 
                            caught[i] = nc;
                            for (t=0; t<*ss; t++)
                                for (l=0; l<*npoly; l++)                            
                                    value[*ss * ((nc-1) * *npoly + l) + t] = 0;
                        }
                        nd++;
                        if (nd >= maxdet) {
                            *resultcode = 2;           /* error */
                            return;  
                        }
                        value[*ss * ((caught[i]-1) * *npoly + k) + s]++;
                        workXY[(nd-1)*2] = xy[0];
                        workXY[(nd-1)*2+1] = xy[1];
                        sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                    }
                }
            }
        }
    }
    for (i=0; i<nd; i++) 
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;    
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingtransect (
    double *lambda,      /* Parameter : expected detection events per metre */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *ntransect,   /* number of different transects */
    int    *kk,          /* number of vertices per transect (vector) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential etc. */
    int    *maxone,      /* maximum of one detection per animal per occasion */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)
{
    int    i,j,k,l,s,t;
    int    nc = 0;
    int    nd = 0;
    int    sumk;
    int    n1,n2;
    int    count;
    int    maxdet;
    int    gotcha;

    int    *cumk;       
    double *cumd;
    struct rpoint *line;
    struct rpoint xy;
    struct rpoint animal;
    double par[3];
    double lx;
    double maxg = 0;
    double lambdak;
    double grx;
    double *workXY;
    int    *sortorder;
    double *sortkey;

    *resultcode = 1; 
    GetRNGstate();

    cumk = (int *) R_alloc(*ntransect+1, sizeof(int));
    cumk[0] = 0;
    for (k =0; k<*ntransect; k++) 
        cumk[k+1] = cumk[k] + kk[k];
    sumk = cumk[*ntransect];
    line = (struct rpoint *) R_alloc(sumk, sizeof(struct rpoint));
    cumd = (double *) R_alloc(sumk, sizeof(double));
    maxdet = *N * *ss * *ntransect * maxperpoly;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));

    /* coordinates of vertices */
    for (i=0; i<sumk; i++) {
        line[i].x = traps[i];
        line[i].y = traps[i+sumk];
    }

    /* cumulative distance along line */
    for (k=0; k<*ntransect; k++) {
        cumd[cumk[k]] = 0;
        for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
            cumd[i+1] = cumd[i] + distance(line[i], line[i+1]);  
        }
    }

    for (i=0; i<*N; i++) caught[i] = 0;                           
    for (s=0; s<*ss; s++) {                       /* each occasion */ 
        for (i=0; i<*N; i++) {                            /* each animal */
            animal.x = animals[i];
            animal.y = animals[i + *N];
            for (k=0; k<*ntransect; k++) {            /* each transect */
                n1 = cumk[k]; 
                n2 = cumk[k+1]-1;
                par[0] = lambda[s];
                par[1] = sigma[s];
                par[2] = z[s];
                lambdak = integral1D (*fn, i, 0, par, 1, traps, animals, n1, n2, sumk, *N);
                if (*maxone)
                    count = (Random() < lambdak);     /* not clear this is useful */
                else
                   count = rpois(lambdak);                          
                maxg = 0;
                if (count>0) {                        /* find maximum - approximate */
                    for (l=0; l<=100; l++) {
                        lx = cumd[n2] * l/100; 
                        xy = getxy (lx, cumd, line, sumk, n1);
                        grx = gr (fn, par, xy, animal);
                        maxg = fmax2(maxg, grx); 
                    }
                    for (l=n1; l<=n2; l++) {
                        xy = line[l];
                        grx = gr (fn, par, xy, animal);
                        maxg = fmax2(maxg, grx); 
                    }
                    maxg= 1.2 * maxg;                 /* safety margin */
                    if (maxg<=0) maxg=0.0001;         /* not found */
                } 
                for (j=0; j<count; j++) {             
                    gotcha = 0;
                    l = 0;
                    while (gotcha == 0) { 
                        lx = Random() * cumd[n2];     /* simulate location */
                        xy = getxy (lx, cumd, line, sumk, n1);
                        grx = gr (fn, par, xy, animal);
                        if (Random() < (grx/maxg))    /* rejection sampling */
                            gotcha = 1;       
                        l++;
                        if (l % 10000 == 0) 
                            R_CheckUserInterrupt();
                        if (l>1e6) gotcha = 1;        /* give up and accept anything!!!! */ 
                    }
                    if (caught[i]==0) {               /* first capture of this animal */
                        nc++; 
                        caught[i] = nc;
                        for (t=0; t<*ss; t++)
                            for (l=0; l<*ntransect; l++)                            
                                value[*ss * ((nc-1) * *ntransect + l) + t] = 0;
                    }
                    nd++;
                    if (nd >= maxdet) {
                        *resultcode = 2;
                        return;  /* error */
                    }
                    value[*ss * ((caught[i]-1) * *ntransect + k) + s]++;
                    workXY[(nd-1)*2] = xy.x;
                    workXY[(nd-1)*2+1] = xy.y;
                    sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                }
            }
        }
    }
    for (i=0; i<nd; i++) 
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;    
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingsignal (
    double *beta0,     /* Parameter : intercept */
    double *beta1,     /* Parameter : slope */
    double *sdS,       /* Parameter : error sd */
    double *cut,       /* detection threshold on transformed scale */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    int    *used,      /* ss x kk array of 0/1 codes for usage */
    int    *fn,        /* code 10 = signal strength, 11 = signal strength with sph spread */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *signal,    /* signal strength, one per detection */  
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)

/* returned signal strength (*fn==10) is on transformed scale */
/* limited to Bernoulli count model binomN = 1 */
{
    double mu, signalvalue;
    int    i,j,k,l,s;
    int    nc = 0;
    int    nd = 0;
    int    maxdet;
    double *worksignal;
    int    *sortorder;
    double *sortkey;

    *resultcode = 1;          
    GetRNGstate();   

    maxdet = *N * *ss * *kk;
    worksignal = (double*) R_alloc(maxdet, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            for (k=0; k<*kk; k++) {
                if (used[s * *kk + k]) {   
                    if (*fn == 10)                     
                        mu  = mufn (i, k, beta0[s], beta1[s], animals, traps, *N, *kk);
                    else
                        mu  = mufnsph (i, k, beta0[s], beta1[s], animals, traps, *N, *kk);
                    signalvalue = norm_rand() * sdS[s] + mu;
                    if (signalvalue > *cut) {                  
                        if (caught[i]==0) {              /* first capture of this animal */
                            nc++; 
                            caught[i] = nc;
                            for (j=0; j<*ss; j++)
                                for (l=0; l<*kk; l++)
                                    value[*ss * ((nc-1) * *kk + l) + j] = 0;
                        } 
                        nd++;
                        if (nd >= maxdet) {
                            *resultcode = 2;
                            return;  /* error */
                        }
                        value[*ss * ((caught[i]-1) * *kk + k) + s] = 1;
                        worksignal[nd-1] = signalvalue;
                        sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                    } 
                }
            }
        }     
    }
    for (i=0; i<nd; i++) 
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) 
        signal[i]   = worksignal[sortorder[i]];
    *n = nc; 
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingtimes (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    int    *used,      /* ss x kk array of 0/1 codes for usage */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *times,     /* time of detection within occasion */  
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)
{
    double lambda, timevalue;
    int    i,j,k,l,s;
    int    nc = 0;
    int    nd = 0;
    double d2val;
    int    maxdet;
    double *work;
    int    *sortorder;
    double *sortkey;

    *resultcode = 1;          
    GetRNGstate();   
    maxdet = *N * *ss * *kk * 100;
    work = (double*) R_alloc(maxdet, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            for (k=0; k<*kk; k++) {
                if (used[s * *kk + k]) {                       
                    timevalue = 0;
                    d2val = d2(i,k, animals, traps, *N, *kk);
                    lambda = pfn(*fn, d2val, g0[s], sigma[s], z[s], 0);
                    if (lambda>0) {
                        while (timevalue < 1) {
                            timevalue += rexp(1/lambda);
                            if (timevalue < 1) {                  
                                if (caught[i]==0) {                 /* first capture of this animal */
                                    nc++; 
                                    caught[i] = nc;
                                    for (j=0; j<*ss; j++)
                                      for (l=0; l<*kk; l++)
                                        value[*ss * ((nc-1) * *kk + l) + j] = 0;
                                } 
                                nd++;
                                if (nd >= maxdet) {
                                    *resultcode = 2;
                                    return;  /* error */
                                }
                                value[*ss * ((caught[i]-1) * *kk + k) + s]++;
                                work[nd-1] = timevalue;
                                sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                            }
                        }
                    }
                }
            }
        }     
    }
    for (i=0; i<nd; i++) 
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) 
        times[i]   = work[sortorder[i]];
    *n = nc; 
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingmulti (
    double *g0,         /* Parameter : detection magnitude  */
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    int    *used,       /* ss x kk array of 0/1 codes for usage */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode  /* 0 for successful completion */
)
  
{
    double *h;
    double hsum[*N];
    double cump[*kk+1];
    double runif;
    int    i,j,k,s;
    int    nc;
    double d2val;
    double p;
  
    *resultcode = 1;
    cump[0] = 0;
    nc = 0;
    GetRNGstate();
    h = (double *) R_alloc(*N * *kk, sizeof(double));
 
    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            hsum[i] = 0;
            for (k=0; k<*kk; k++)
            {
                d2val = d2(i,k, animals, traps, *N, *kk);
                p = pfn(*fn, d2val, g0[s], sigma[s], z[s], 0);  
                p = p * used[s * *kk + k];           /* zero if not used 2009 11 09 */
                h[k * *N + i] = -log(1 - p);
                hsum[i] += h[k * *N + i];
            }

            for (k=0; k<*kk; k++) {
                cump[k+1] = cump[k] + h[k * *N + i]/hsum[i];
            }

            if (Random() < (1-exp(-hsum[i])))
            {
               if (caught[i]==0)           /* first capture of this animal */
               {
                   nc++; 
                   caught[i] = nc;
                   for (j=0; j<*ss; j++) 
                       value[*ss * (nc-1) + j] = 0;
               } 
               runif = Random();
               k = 0;
               while ((runif > cump[k]) & (k<*kk)) k++;  /* pick a trap */
               value[*ss * (caught[i]-1) + s] = k;
            }
        } 
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();

}
/*==============================================================================*/

void trappingsingle (
    double *g0,        /* Parameter : detection magnitude  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    int    *used,      /* ss x kk array of 0/1 codes for usage */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode /* 0 for successful completion */
)
{
    int    i,j,k,s;
    int    nc         = 0;
    int    tr_an_indx = 0;
    double d2val;
    double p;
    int nanimals;         /* temporary */
    int ntraps;           /* temporary */
    int occupied[*kk];    /* today */
    int intrap[*N];       /* today   */
  
    struct  trap_animal *tran;  
    double event_time;
    int anum;
    int tnum;
    int nextcombo;
    int finished;
    int OK;

    /* MAIN LINE */
    *resultcode = 1;
    GetRNGstate();
    tran = (struct trap_animal *) R_alloc(*N * *kk, sizeof(struct trap_animal));
    for (i=0; i<*N; i++) caught[i] = 0;   /* has animal i been caught in session? */
    for (s=0; s<*ss; s++) {
        /* initialise day */
        tr_an_indx = 0;
        nanimals = *N;
        ntraps   = *kk;
        for (i=0; i<*N; i++) intrap[i] = 0;
        for (k=0; k<*kk; k++) occupied[k] = 0;
        nextcombo = 0;

        /* make tran */
        for (i=0; i<*N; i++)   /* animals */
        for (k=0; k<*kk; k++)  /* traps */
        if (used[s * *kk + k]) {
            d2val = d2(i,k, animals, traps, *N, *kk);
            p = pfn(*fn, d2val, g0[s], sigma[s], z[s], 0);         
            event_time = randomtime(p);
            if (event_time <= 1) {
                tran[tr_an_indx].time   = event_time;
                tran[tr_an_indx].animal = i;    /* 0..*N-1 */
                tran[tr_an_indx].trap   = k;    /* 0..*kk-1 */
                tr_an_indx++;
            }
        }

        if (tr_an_indx>0) probsort (tr_an_indx, tran);

        /* make captures */
        while ((nextcombo < tr_an_indx) & (nanimals>0) & (ntraps>0)) {
            finished = 0;
            OK       = 0; 
            while ((1-finished)*(1-OK) > 0) {    /* until finished or OK */
                if (nextcombo >= (tr_an_indx)) finished = 1;  /* no more to process */
                else {
                    anum = tran[nextcombo].animal;
                    tnum = tran[nextcombo].trap;
                    OK = (1-occupied[tnum]) * (1-intrap[anum]); /* not occupied and not intrap */
                    nextcombo++;       
                }
            }
            if (finished==0) {                   /* Record this capture */
                  occupied[tnum] = 1;
                  intrap[anum]   = tnum+1;       /* trap = k+1 */
                  nanimals--;
                  ntraps--;
            }
        }
        for (i=0; i<*N; i++) 
        if (intrap[i]>0) {
            if (caught[i]==0) {                  /* first capture of this animal */
               nc++; 
               caught[i] = nc;                   /* nc-th animal to be captured */
               for (j=0; j<*ss; j++) 
                   value[*ss * (nc-1) + j] = 0;
             } 
             value[*ss * (caught[i]-1) + s] = intrap[i];  /* trap = k+1 */
        } 
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();

}
/*==============================================================================*/

int rdiscrete (int n, double pmix[]) 
/* return random discrete observation from distribution in pmix */
{
    double *cumpmix;
    int x;
    double r;
    if (n<1) error ("invalid n in rdiscrete");
    if (n==1) return (0);
    else {
        cumpmix = (double*) R_alloc(n+1, sizeof(double));   
        cumpmix[0] = 0;
        for (x=0; x<n; x++) {
            cumpmix[x+1] = cumpmix[x] + pmix[x];
        }
        r = Random();
        for (x=1; x<=n; x++) if (r<cumpmix[x]) break;
        return(x);
    }
}

void simsecr (
    int    *detect,     /* detector 0 multi, 1 proximity, 2 single, 3 count, 4 area ??? */ 
    double *gsb0val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    double *gsb1val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [caught before] */
    int    *cc0,        /* number of g0/sigma/b combinations for naive animals */
    int    *cc1,        /* number of g0/sigma/b combinations for caught before */
    int    *gsb0,       /* lookup which g0/sigma/b combination to use for given g, S, K [naive animal] */
    int    *gsb1,       /* lookup which g0/sigma/b combination to use for given n, S, K [caught before] */
    int    *N,          /* number of animals */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *nmix,       /* number of classes */
    double *animals,    /* x,y points of animal range centres (first x, then y) */
    double *traps,      /* x,y locations of traps (first x, then y) */
    int    *used,       /* ss x kk array of 0/1 codes for usage */
    int    *Markov,     /* code 0 if behavioural response is learned, 1 if Markov */
    int    *binomN,     /* number of trials for 'count' detector modelled with binomial */
    double *cut,        /* detection threshold on transformed scale */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 = uniform */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* sequence number in session (0 if not caught) */
    double *detectedXY, /* x,y locations of detections  */
    double *signal,     /* vector of signal strengths, one per detection */
    int    *value,      /* return value array of trap locations n x s */
    int    *resultcode
)
{

    double d2val;
    double p;
    int    i,j,k,l,s;
    int    nc = 0;
    int    nk = 0;             /* number of detectors (polygons or transects when *detect==6,7) */
    int    count = 0;
    int    *caughtbefore;
    int    *x;                 /* mixture class of animal i */
    double *pmix;
    double runif;
    int    wxi = 0;
    int    c = 0;
    int    gpar = 2;
    double g0 = 0;
    double sigma = 0;
    double z = 0;
    double *work = NULL;
    int    *sortorder = NULL;
    double *sortkey = NULL;

    /*
        *detect may take values -
       -1  single-catch traps
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  binary quadrat search
        4  count  quadrat search
        5  signal detectors
        6  polygon detectors 
        7  transect detectors 
    */
    /*========================================================*/
    /* 'single-catch only' declarations */
    int    tr_an_indx = 0;  
    int    nanimals;        
    int    ntraps;          
    int    *occupied = NULL;   
    int    *intrap = NULL;      
    struct trap_animal *tran = NULL;  
    double event_time;
    int    anum;
    int    tnum;
    int    nextcombo;
    int    finished;
    int    OK;

    /*========================================================*/
    /* 'multi-catch only' declarations */
    double *h = NULL;          /* multi-catch only */
    double *hsum = NULL;       /* multi-catch only */
    double *cump = NULL;       /* multi-catch only */

    /*========================================================*/
    /* 'polygon & transect only' declarations */
    int    nd = 0; 
    int    cumk[maxnpoly+1]; 
    int    sumk;               /* total number of vertices */
    int    g=0;
    int    *gotcha;
    double xy[2]; 
    int    n1,n2,t;
    double par[3];
    int    np = 1;             /* n points each call of gxy */
    double w;
    int    maxdet=1;
    double *cumd = NULL;
    struct rpoint *line = NULL;
    struct rpoint xyp;
    struct rpoint animal;
    double lx;
    double maxg = 0;
    double lambdak;  /* temp value for Poisson rate */
    double grx;      /* temp value for integral gr */

    /*========================================================*/
    /* 'signal-strength only' declarations */
    double beta0;
    double beta1;
    double sdS;
    double mu;
    double signalvalue;

    /*========================================================*/
    /* MAIN LINE */

    gotcha = &g;
    *resultcode = 1;
    caughtbefore = (int *) R_alloc(*N, sizeof(int));
    x = (int *) R_alloc(*N, sizeof(int));
    for (i=0; i<*N; i++) x[i] = 0;
    pmix = (double *) R_alloc(*nmix, sizeof(double));
    if ((*detect < -1) | (*detect > 7)) return;
    if (*detect == -1) {                                   /* single-catch only */
        occupied = (int*) R_alloc(*kk, sizeof(int));
        intrap = (int*) R_alloc(*N, sizeof(int));
        tran = (struct trap_animal *) R_alloc(*N * *kk,  
            sizeof(struct trap_animal));
    }
    if (*detect == 0) {                                    /* multi-catch only */    
        h = (double *) R_alloc(*N * *kk, sizeof(double));  
        hsum = (double *) R_alloc(*N, sizeof(double));
        cump = (double *) R_alloc(*kk+1, sizeof(double));
        cump[0] = 0;                                     
    }
    if (*detect == 5) {                                    /* signal only */
        maxdet = *N * *ss * *kk;                           
        if (!((*fn == 10) | (*fn == 11)))
            error ("simsecr not implemented for this combination of detector & detectfn");

    }
    if ((*detect == 6) | (*detect == 7)) {                 /* polygon or transect only */
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {                       /* maxnpoly much larger than npoly */ 
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
        sumk = cumk[nk];
        maxdet = *N * *ss * nk * maxperpoly;
    }
    else nk = *kk;
    if (*detect == 7) {                                    /* transect only */
        line = (struct rpoint *) R_alloc(sumk, sizeof(struct rpoint));
        cumd = (double *) R_alloc(sumk, sizeof(double));
        /* coordinates of vertices */
        for (i=0; i<sumk; i++) {
            line[i].x = traps[i];
            line[i].y = traps[i+sumk];
        }
        /* cumulative distance along line; all transects end on end */
        for (k=0; k<nk; k++) {
            cumd[cumk[k]] = 0;
            for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
                cumd[i+1] = cumd[i] + distance(line[i], line[i+1]);  
            }
        }
    }    
    if ((*detect==5) | (*detect==6) | (*detect==7)) {
        work = (double*) R_alloc(maxdet*2, sizeof(double));   /* twice size needed for signal */
        sortorder = (int*) R_alloc(maxdet, sizeof(int));
        sortkey = (double*) R_alloc(maxdet, sizeof(double));
    }
    GetRNGstate();

    /* may be better to pass pmix */
    gpar = 2;
    if ((*fn == 1) | (*fn == 3) | (*fn == 5) | (*fn == 10) | (*fn == 11)) gpar ++;
    if (*nmix>1) gpar++;



    if (*nmix>1) {
        if (*nmix>2) 
            error("simsecr nmix>2 not implemented");
        for (i=0; i<*nmix; i++) {
            wxi = i4(0,0,0,i,*N,*ss,nk);
            c = gsb0[wxi] - 1; 
            pmix[i] = gsb0val[*cc0 * (gpar-1) + c];    /* assuming 4-column gsb */
        }
        for (i=0; i<*N; i++) {
            x[i] = rdiscrete(*nmix, pmix) - 1;
        }
    }

    for (i=0; i<*N; i++) {
        caught[i] = 0;
    }

    /* MAIN LOOP */
    for (s=0; s<*ss; s++) {

        /* --------------------------------------------- */
        /* universal update of 'previous-capture' status */
        for (i=0; i<*N; i++) {
            if ((s>0) & *Markov) {
                caughtbefore[i] = 0;
                for (k=0;k<nk;k++)
                    if (value[i3(s-1, k, i, *ss, nk)]) 
                        caughtbefore[i] = caught[i];
            }
            else
                caughtbefore[i] = caught[i];
        }

        /* ------------------ */
        /* single-catch traps */
        if (*detect == -1) {  
            /* initialise day */
            tr_an_indx = 0;
            nanimals = *N;
            ntraps   = nk;
            for (i=0; i<*N; i++) intrap[i] = 0;
            for (k=0; k<nk; k++) occupied[k] = 0;
            nextcombo = 0;

            /* make tran */
            for (i=0; i<*N; i++) {  /* animals */
                for (k=0; k<nk; k++) { /* traps */
                    if (used[s * nk + k]) {
                        wxi = i4(i,s,k,x[i],*N,*ss,nk);
                        if (caughtbefore[i]==0) {
                            c = gsb0[wxi]-1;
                            g0 = gsb0val[c];
                            sigma = gsb0val[*cc0 + c];
                            if ((*fn==1) | (*fn == 5)) z = gsb0val[2* *cc0 + c];
                        }
                        else {                 
                            c = gsb1[wxi]-1;
                            g0 = gsb1val[c];
                            sigma = gsb1val[*cc1 + c];
                            if ((*fn==1) | (*fn == 5)) z = gsb1val[2* *cc1 + c];
                        }                       
                        d2val = d2(i,k, animals, traps, *N, nk);
                        p = pfn(*fn, d2val, g0, sigma, z, 0);         
                        event_time = randomtime(p);
                        if (event_time <= 1) {
                            tran[tr_an_indx].time   = event_time;
                            tran[tr_an_indx].animal = i;    /* 0..*N-1 */
                            tran[tr_an_indx].trap   = k;    /* 0..nk-1 */
                            tr_an_indx++;
                        }
                    }
                }
            }
            /* end of make tran */

            if (tr_an_indx > 1) probsort (tr_an_indx, tran);        

            while ((nextcombo < tr_an_indx) & (nanimals>0) & (ntraps>0)) {
                    finished = 0;
                OK       = 0; 
                    while ((1-finished)*(1-OK) > 0) {      /* until finished or OK */
                    if (nextcombo >= (tr_an_indx)) 
                            finished = 1;                  /* no more to process */
                    else {
                            anum = tran[nextcombo].animal;
                        tnum = tran[nextcombo].trap;
                            OK = (1-occupied[tnum]) * (1-intrap[anum]); /* not occupied and not intrap */
                        nextcombo++;       
                        }
                }             
                    if (finished==0) {
                       /* Record this capture */
                          occupied[tnum] = 1;
                      intrap[anum]   = tnum+1;         /* trap = k+1 */
                          nanimals--;
                      ntraps--;
                    }
            }    
                for (i=0; i<*N; i++) {
                if (intrap[i]>0) {
                    if (caught[i]==0) {                    /* first capture of this animal */
                       nc++; 
                       caught[i] = nc;                     /* nc-th animal to be captured */
                       for (j=0; j<*ss; j++) 
                           value[*ss * (nc-1) + j] = 0;
                    } 
                    value[*ss * (caught[i]-1) + s] = intrap[i];  /* trap = k+1 */
                } 
            }
        }

        /* -------------------------------------------------------------------------- */
        /* multi-catch trap; only one site per occasion (drop last dimension of capt) */
        else if (*detect == 0) {  
            for (i=0; i<*N; i++) {
                hsum[i] = 0;
                for (k=0; k<nk; k++) {
                    wxi = i4(i,s,k,x[i],*N,*ss,nk);
                    if (caughtbefore[i]==0) {
                        c = gsb0[wxi]-1;
                        g0 = gsb0val[c];
                        sigma = gsb0val[*cc0 + c];
                        if ((*fn==1) | (*fn == 5)) z = gsb0val[2* *cc0 + c];
                    }
                    else {                 
                        c = gsb1[wxi]-1;
                        g0 = gsb1val[c];
                        sigma = gsb1val[*cc1 + c];
                        if ((*fn==1) | (*fn == 5)) z = gsb1val[2* *cc1 + c];
                    }
                    d2val = d2(i,k, animals, traps, *N, nk);
                    p = pfn(*fn, d2val, g0, sigma, z, 0);         
                    p = p * used[s * nk + k];           /* zero if not used */
                    h[k * *N + i] = -log(1 - p);
                    hsum[i] += h[k * *N + i];
                }

                for (k=0; k<nk; k++) {
                    cump[k+1] = cump[k] + h[k * *N + i]/hsum[i];
                }                 
                if (Random() < (1-exp(-hsum[i]))) {
                    if (caught[i]==0)  {        /* first capture of this animal */
                        nc++; 
                        caught[i] = nc;
                        for (j=0; j<*ss; j++)
                            value[*ss * (nc-1) + j] = 0;
                    } 
                    /* find trap with probability proportional to p
                       searches cumulative distribution of p  */
                    runif = Random();
                    k = 0;
                    while ((runif > cump[k]) & (k<nk)) k++;
                    value[*ss * (caught[i]-1) + s] = k;  /* trap = k+1 */
                }
            }
        }

        /* -------------------------------------------------------------------------------- */
        /* the 'proximity' group of detectors 1:4 - proximity, count, quadratbinary, quadratcount */
        else if ((*detect >= 1) & (*detect <= 4)) {
            for (i=0; i<*N; i++) {
                for (k=0; k<nk; k++) {
                    if (used[s * nk + k]) {
                        wxi = i4(i,s,k,x[i],*N,*ss,nk);
                        if (caughtbefore[i]==0) {
                            c = gsb0[wxi]-1;
                            g0 = gsb0val[c];
                            sigma = gsb0val[*cc0 + c];
                            if ((*fn==1) | (*fn == 5)) z = gsb0val[2* *cc0 + c];
                        }
                        else {                 
                            c = gsb1[wxi]-1;
                            g0 = gsb1val[c];
                            sigma = gsb1val[*cc1 + c];
                            if ((*fn==1) | (*fn == 5)) z = gsb1val[2* *cc1 + c];
                        }
          
                        d2val = d2(i,k, animals, traps, *N, nk);
                        p = pfn(*fn, d2val, g0, sigma, z, 0);         
    
                        if (p < -0.1) { PutRNGstate(); return; }   /* error */   
    
                        if (p>0) {
                            if (*detect == 1) {
                                count = Random() < p;              /* binary proximity */
                            }
                            else if (*detect == 2) {               /* count proximity */
                                if (*binomN == 0)
                                    count = rpois(p);
                                else if (*binomN == 1)
                                    count = Random() < p;              
                                else if (*binomN < 0) {
                                    /* prob = size / (size + mu) */
                                    *binomN = abs(*binomN);
                                    count = rnbinom(*binomN, *binomN / (*binomN+p));   
                                }
                                else
                                    count = rbinom(*binomN, p / *binomN);   
                            }
                            else if (*detect == 3) {               /* quadrat binary */
                                count = Random() < 1-(exp(-p));    /* Poisson Pr(nonzero) */
                            }
                            else if (*detect == 4) {               /* quadrat count */
                                count = rpois(p);
                            }         
                            if (count>0) {
                                if (caught[i]==0) {                /* first capture of this animal */
                                    nc++; 
                                    caught[i] = nc;
                                    for (j=0; j<*ss; j++)
                                      for (l=0; l<nk; l++)
                                        value[*ss * ((nc-1) * nk + l) + j] = 0;
                                } 
                                value[*ss * ((caught[i]-1) * nk + k) + s] = count;
                            }
                        }
                    }
                }
            } 
        }    

        /* -------------------------------------------------------------------------------- */
        /* polygon detectors  */
        else if (*detect == 6) {
            for (i=0; i<*N; i++) {
                /* this implementation assumes NO VARIATION AMONG DETECTORS */
                wxi = i4(i,s,0,x[i],*N,*ss,nk);
                if (caughtbefore[i]==0) {
                    c = gsb0[wxi]-1;
                    g0 = gsb0val[c];
                    sigma = gsb0val[*cc0 + c];
                    if ((*fn==1) | (*fn == 5)) 
                        z = gsb0val[2* *cc0 + c];
                }
                else {                 
                    c = gsb1[wxi]-1;
                    g0 = gsb1val[c];
                    sigma = gsb1val[*cc1 + c];
                    if ((*fn==1) | (*fn == 5)) 
                        z = gsb1val[2* *cc1 + c];
                }
                count = rpois(g0);                              /* number of cues */
                w = 5 * sigma;
                par[0] = 1;
                par[1] = sigma;
                par[2] = z;
                for (j=0; j<count; j++) {             
                    gxy (&np, fn, par, &w, xy);                 /* simulate location */
                    xy[0] = xy[0] + animals[i];
                    xy[1] = xy[1] + animals[*N + i];
                    for (k=0; k<nk; k++) {                      /* each polygon */
                        if (used[s * nk + k]) {
                            n1 = cumk[k]; 
                            n2 = cumk[k+1]-1;
                            inside(xy, &n1, &n2, &sumk, traps, gotcha);  /* assume closed */
                            if (*gotcha > 0) {
                                if (caught[i]==0) {             /* first capture of this animal */
                                    nc++; 
                                    caught[i] = nc;
                                    for (t=0; t<*ss; t++)
                                        for (l=0; l<nk; l++)                            
                                            value[*ss * ((nc-1) * nk + l) + t] = 0;
                                }
                                nd++;
                                if (nd > maxdet) {
                                    *resultcode = 2;
                                    return;  /* error */
                                }
                                value[*ss * ((caught[i]-1) * nk + k) + s]++;
                                work[(nd-1)*2] = xy[0];
                                work[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                            }
                        }
                    }
                }
            }
        }
    
        /* -------------------------------------------------------------------------------- */
        /* transect detectors  */
        else if (*detect == 7) {
            for (i=0; i<*N; i++) {                            /* each animal */
                animal.x = animals[i];
                animal.y = animals[i + *N];
                for (k=0; k<nk; k++) {                        /* each transect */
                    if (used[s * nk + k]) {
                        wxi = i4(i,s,k,x[i],*N,*ss,nk);
                        if (caughtbefore[i]==0) {
                            c = gsb0[wxi]-1;
                            par[0] = gsb0val[c];              /* lambda(0) */
                            par[1] = gsb0val[*cc0 + c];       /* sigma */
                            if ((*fn==1) | (*fn == 5))        /* z */
                                par[2] = gsb0val[2* *cc0 + c];
                        }
                        else {                 
                            c = gsb1[wxi]-1;
                            par[0] = gsb1val[c];              /* naive lambda(0) */
                            par[1] = gsb1val[*cc1 + c];       /* naive sigma */
                            if ((*fn==1) | (*fn == 5))        /* naive z */
                                par[2] = gsb1val[2* *cc1 + c];
                        }
                        n1 = cumk[k]; 
                        n2 = cumk[k+1]-1;
                        lambdak = integral1D (*fn, i, 0, par, 1, traps, animals, n1, n2, sumk, *N);
                        count = rpois(lambdak);               /* number of detections on this transect */               
                        maxg = 0;
                        if (count>0) {                        /* find maximum - approximate */
                            for (l=0; l<=100; l++) {
                                lx = cumd[n2] * l/100; 
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = fmax2(maxg, grx); 
                            }
                            for (l=n1; l<=n2; l++) {
                                xyp = line[l];
                                grx = gr (fn, par, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = fmax2(maxg, grx); 
                            }
                            maxg= 1.2 * maxg;                 /* safety margin */
                            if (maxg<=0) 
                                Rprintf("maxg error in simsecr\n"); /* not found */

                        } 
                        for (j=0; j<count; j++) {             
                            *gotcha = 0;
                            l = 0;
                            while (*gotcha == 0) { 
                                lx = Random() * cumd[n2];     /* simulate location */
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xyp, animal);
                                if (Random() < (grx/maxg))   /* rejection sampling */
                                    *gotcha = 1;       
                                l++;
                                if (l % 10000 == 0) 
                                    R_CheckUserInterrupt();
                            }
                            if (caught[i]==0) {               /* first capture of this animal */
                                nc++; 
                                caught[i] = nc;
                                for (t=0; t<*ss; t++)
                                    for (l=0; l<nk; l++)                            
                                        value[*ss * ((nc-1) * nk + l) + t] = 0;
                            }
                            nd++;
                            if (nd >= maxdet) {
                                *resultcode = 2;  /* error */
                                return;  
                            }
                            value[*ss * ((caught[i]-1) * nk + k) + s]++;
                            work[(nd-1)*2] = xyp.x;
                            work[(nd-1)*2+1] = xyp.y;
                            sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                        }
                    } 
                }                                             /* end loop over transects */
            }                                                 /* end loop over animals */
        }
   
        /* ------------------------ */
        /* signal strength detector */
        else if (*detect == 5) {
            for (i=0; i<*N; i++) {
                for (k=0; k<nk; k++) {
                    if (used[s * nk + k]) {
                        wxi = i4(i,s,k,x[i],*N,*ss,nk);
                        c = gsb0[wxi]-1;
                        beta0 = gsb0val[c];
                        beta1 = gsb0val[*cc0 + c];
                        sdS   = gsb0val[2* *cc0 + c];
                        if (*fn == 10)
                            mu  = mufn (i, k, beta0, beta1, animals, traps, *N, nk);
                        else 
                            mu  = mufnsph (i, k, beta0, beta1, animals, traps, *N, nk);
                        signalvalue = norm_rand() * sdS + mu;
                        if (signalvalue > *cut) {
                            if (caught[i]==0) {                /* first capture of this animal */
                                nc++; 
                                caught[i] = nc;
                                for (j=0; j<*ss; j++)
                                  for (l=0; l<nk; l++)
                                    value[*ss * ((nc-1) * *kk + l) + j] = 0;
                            } 
                            nd++;
                            value[*ss * ((caught[i]-1) * *kk + k) + s] = 1;
                            work[nd-1] = signalvalue;
                            sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                        }
                    }
                }
            }
        } 
    }   /* loop over s */

    if ((*detect==5) | (*detect==6) | (*detect==7)) {
        for (i=0; i<nd; i++) sortorder[i] = i;
        if (nd>0) rsort_with_index (sortkey, sortorder, nd);
        if (*detect==5) 
            for (i=0; i<nd; i++) signal[i] = work[sortorder[i]];
        else {
            for (i=0; i<nd; i++) {
                detectedXY[i]    = work[sortorder[i]*2];
                detectedXY[i+nd] = work[sortorder[i]*2+1];
            }
        }
    }
    *n = nc;
    PutRNGstate();
    *resultcode = 0;
}

/*==============================================================================*/

/*
    'naive' functions are used to estimate auto initial values
    these use only the halfnormal detection function 4/5/08
*/

void naived (
  double *sigma,   /* Parameter : detection scale */
  int    *kk,      /* number of traps */
  int    *nc,
  double *traps,   /* x,y locations of traps (first x, then y)   */
  double *animals, /* x,y locations of animals (first x, then y) */
  int    *fn,      /* code 0 = halfnormal ONLY */
  double *value    /* return value */
)
{
  double truncate2 = (2.45 * *sigma) * (2.45 * *sigma);
  double sump  = 0;
  double sumdp = 0;
  double x,y;
  double dij, d21, d22, p1p2;
  int i,j,n;

  for (n=0; n<*nc; n++)
  {
    x = animals[n];
    y = animals[n + *nc];

    for (i=0; i<*kk; i++)
      for (j=0; j<(i-1); j++)
        {

        dij = (traps[i] - traps[j]) * (traps[i] - traps[j]) + 
                (traps[i+*kk] - traps[j+*kk]) * (traps[i+*kk] - traps[j+*kk]);
        d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+*kk] - y) * (traps[i+*kk] - y);
        d22 = (traps[j] - x) * (traps[j] - x) + (traps[j+*kk] - y) * (traps[j+*kk] - y);

        if ((d21<=truncate2) & (d22<=truncate2)) 
           p1p2 = exp(-(d21+d22) / 2 / *sigma / *sigma);
        else
           p1p2 = 0;

        sump  += p1p2;
        sumdp += p1p2 * sqrt(dij);
  
        }
    for (i=0; i<*kk; i++)  /* diagonal */
      {
        d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+*kk] - y) * (traps[i+*kk] - y);
        if (d21<=truncate2)                                     /* d21=d22 */
          sump += exp(-2*d21 /2 / *sigma / *sigma)/2;
      }
  }
  *value = sumdp/sump;
}
/*==============================================================================*/

void naiveRPSV (
  double *sigma,   /* Parameter : detection scale */
  int    *kk,      /* number of traps */
  int    *nc,
  double *traps,   /* x,y locations of traps (first x, then y)   */
  double *animals, /* x,y locations of animals (first x, then y) */
  int    *fn,      /* code 0 = halfnormal ONLY */
  double *value    /* return value */
)
{
    double truncate2 = (2.45 * *sigma) * (2.45 * *sigma);
    double x,y;
    int k, n;
    double d2k;
    double pk;
    double pdot;
    double sumd2k;
    double sumpk;  
    double sump  = 0;
    double sumpRPSV = 0;
    for (n=0; n<*nc; n++)
    {
        x = animals[n];
        y = animals[n + *nc];  
        pdot   = 1;
        sumd2k = 0;
        sumpk  = 0;
        for (k=0; k< *kk; k++) {
            d2k = (traps[k] - x) * (traps[k] - x) + (traps[k + *kk]-y) * (traps[k + *kk]-y);
            if (d2k <= truncate2)
                pk = exp(-d2k / 2 / *sigma / *sigma);
            else 
                pk = 0;
            sumd2k += pk * d2k;
            sumpk  += pk;
            pdot   *= (1 - pk);
        }
        pdot     = 1 - pdot;                  /* overall detection probability (excl g0) */
        sump     += pdot;
    
        if (sumd2k/sumpk > 0)
            sumpRPSV += pdot * sqrt(sumd2k/sumpk);
    }

    if (sump > fuzz)
        *value = sumpRPSV/sump;
    else 
        *value = -1;

}

/*==============================================================================*/

void naivecap2 (
  int    *detect,  /* code 0 = multicatch, 1 = proximity */
  double *g0,      /* Parameter : detection magnitude */
  double *sigma,   /* Parameter : detection scale */
  int    *ss,      /* number of occasions */
  int    *kk,      /* number of traps */
  int    *mm,
  double *traps,   /* x,y locations of traps (first x, then y)  */
  double *mask,    /* x,y points on mask (first x, then y)  */
  int    *fn,      /* code 0 = halfnormal ONLY */
  double *value    /* return value  */
)
{
    double product;
    double d2val;
    double pk;
    int m,k;
    double nsum = 0;
    double psum = 0;
  
    for (m=0; m<*mm; m++)
    {
        product = 1.0;
        for (k=0; k<*kk; k++)
        {
            d2val = d2(m, k, mask, traps, *mm, *kk);
            pk = *g0 * expmin(-d2val / 2 / *sigma / *sigma);
            product *= (1 - pk);
            if (*detect == 1) nsum += pk;
        }      
        if (*detect == 0) nsum += (1 - product);
        psum += 1 - pow(product, *ss);
    }
    if (psum<=0)
      *value = 0;    /* failed */
    else
      *value = *ss * nsum / psum;
}
/*==============================================================================*/

void makelookup (   
  double *x,            /* input matrix */
  int    *nrow,         /* input */
  int    *ncol,         /* input */
  int    *unique,       /* output number of unique rows */
  double *y,            /* output matrix of unique rows (byrow=T) */
  int    *index,        /* output lookup rows of x in y */
  int    *resultcode)   /* zero if OK */

/*
   Create lookup table to the unique rows in a matrix
   Only the first 'unique' rows of y contain valid data on exit
   The indices are 1..length(unique(x))=nrow(y')
   MGE 2008 05 07
*/

{
    int i;    
    int j;
    int k;
    int dupl = 0;
    double *ytemp;

    *resultcode = 1;

    ytemp = (double *) R_alloc(*nrow * *ncol, sizeof (double));

    /*
        avoid sort for now as it's complex to keep order of first occurrence, not needed
        scan for unique rows of x, copying to y
        assign unique index to original rows as we go    
    */

    for (j=0; j < *ncol; j++) ytemp[j] = x[*nrow * j];    /* first row */
    index[0] = 1;
    *unique=0;

    for (i=1; i < *nrow; i++) {
       /* Is this row unique? Compare with each previous unique in turn */
       for (k=0; k <= *unique; k++) {
           dupl = 1;
           for (j=0; j < *ncol; j++) {
               if (x[*nrow * j + i] != ytemp[k * *ncol + j]) 
               {dupl=0; break;}
           }
           if (dupl==1) break;  /* found previous instance */
       }
       if (dupl==0) { /* add unique row */
           *unique = *unique + 1;
           k = *unique;
           for (j=0; j< *ncol; j++) ytemp[k * *ncol + j] = x[*nrow * j + i];          
       }
       index[i] = k+1;
    }

    *unique = *unique + 1;   /* number of unique rows */
    for (i=0; i<(*unique * *ncol); i++) y[i] = ytemp[i];

    *resultcode = 0;

}
/*==============================================================================*/

void nearest (   
  double *xy,       /* input point */
  int    *ntrap,    /* input */
  double *traps,    /* input */
  int    *p,        /* output index of nearest point */
  double *d)        /* output distance to nearest point */
{
    int i;
    int id=-1;
    double d2;
    *d = 1e100;
    for (i=0; i<*ntrap; i++) 
    {
        d2 = (traps[i] - xy[0]) * (traps[i] - xy[0]) +
             (traps[i + *ntrap] - xy[1]) * (traps[i + *ntrap] - xy[1]);
        if (d2 < *d) { *d = d2; id = i; }
    }
    *d = sqrt(*d);  
    *p = id+1;
}
/*==============================================================================*/

void inside (
    double *xy, 
    int    *n1, 
    int    *n2, 
    int    *npts,
    double *poly, 
    int    *in)
{
/*
    Is point xy inside poly?
    Based on contribution on s-news list by Peter Perkins 23/7/96
    We assume poly is closed, and in col-major order (x's then y's)
*/
  
    double theta = 0;
    double cutoff = 1e-6;
    int k;
    int ns;
    double N;
    double d;
    double *temp;

    ns = *n2 - *n1 + 1;   /* number of selected points */
    temp = (double *) R_alloc((ns+1) * 2, sizeof (double));
    
    /* get & translate to coords centered at each test point */
    for (k=0; k < ns; k++) 
    {
        temp[k]     = poly[k + *n1] - xy[0];           /* x */
        temp[k + ns] = poly[k + *n1 + *npts] - xy[1];    /* y */
    }
  
    for (k=0; k < (ns-1); k++) 
    {
        N = temp[k] * temp[k+1 + ns] - temp[k + ns] * temp[k+1];
        d = temp[k] * temp[k+1]      + temp[k + ns] * temp[k+1 + ns];  
        if (abs(d)>0) { N = N/abs(d);  d = d/abs(d); }
        theta += atan2(N, d);
    }
    theta = abs(theta);
    if (abs(theta - 2* M_PI) < cutoff)    /* M_PI is Rmath.h constant */
        *in = 1;    /* inside */
    else
        *in = 0;    /* outside */
}
/*==============================================================================*/

void ontransect (
    double *xy, 
    int    *n1, 
    int    *n2, 
    int    *npts,
    double *transect, 
    double *tol,
    int    *on)
{
/*
    Is point xy on transect?
    We assume transect coordinates are in col-major order (x's then y's)
    http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/ 29/11/09
*/ 
    int k;
    double r;
    double u;
    struct rpoint p,p1,p2,p3;
    double minr = 1e20;
    *on = 0;

    p3.x = xy[0];
    p3.y = xy[1];

    for (k= *n1; k < *n2; k++) 
    {        
        p1.x = transect[k];
        p1.y = transect[k+*npts];
        p2.x = transect[k+1];
        p2.y = transect[k+1+*npts];
        if (distance(p1,p2) > 0) { 
            u = ((p3.x-p1.x) * (p2.x-p1.x) + (p3.y-p1.y) * (p2.y-p1.y)) / 
                ((p2.x-p1.x) * (p2.x-p1.x) + (p2.y-p1.y) * (p2.y-p1.y));
            if ((u>=0) && (u<=1)) {
                p.x = p1.x + u * (p2.x-p1.x);
                p.y = p1.y + u * (p2.y-p1.y);
                r = distance (p,p3);
                minr = fmin2(r,minr);
            }
        }
    }
    /* check at each vertex to be sure */
    for (k= *n1; k <= *n2; k++) {
        p1.x = transect[k];
        p1.y = transect[k+*npts];
        r = distance (p1,p3);
        minr = fmin2(r,minr);
    }
    if (minr < *tol)
        *on = 1;    /* on transect */
    else
        *on = 0;    /* off transect */
}
/*==============================================================================*/
