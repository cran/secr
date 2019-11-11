#include "secr.h"
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>

// C code for R package 'secr' 

// Murray Efford 

// This file contains functions for integration of home range overlap
// with polygon detectors and similar.
// 
// The choice of detection function is a bit problematic. Before Jan 2016 (2.10.1) 
// the function integrated was g(r) (a probability, hence between 0 and 1) but the theory 
// is for the hazard. Testing zfn 2016-01-02.
// 
// Shifting entirely to hazard formulation (zfn) 2017-03-22

// find upper and lower points on perimeter of poly at x-coordinate x 
void yab(double x[], int *i, int *np, double poly[], double *a, double *b) {
    int k;
    int nv = 0;
    double ab[3];
    // note 'sign' is RMath function 
    for (k=0; k< (*np-1); k++) {
        if (R::sign(poly[k]- x[*i]) != R::sign(poly[k+1]- x[*i])) {
           ab[nv] = poly[k+ *np] + (x[*i]-poly[k]) * (poly[k+1+*np]-poly[k+*np]) /
              (poly[k+1]-poly[k]);
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
    fnptr fnzr; // = zhnr;
    NumericVector gsb(3);   // passed to fnzr
    p = (double*) ex;
    gsb[0] = p[0];
    gsb[1] = p[1];
    gsb[2] = p[2];
    fn = round(p[3]);
    mx = p[4];
    my = p[5];
    xy[0] = p[6];
    // set detection function 
    fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
    for (i=0; i<n; i++) {
        xy[1] = x[i];   // set each y value 
        d = sqrt ( (xy[1]-my)*(xy[1]-my) + (xy[0]-mx)*(xy[0]-mx) );
        x[i] = fnzr(gsb, d);   // z(r) 
    }
}

void fx(double *x, int n, void *ex) {
    int i;
    double * p;

    double poly[maxvertices * 2];
    double a;
    double b;

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
    int kk;
    p = (double*) ex;
    kk = round(p[9]);

    for (i=0; i<kk; i++) {
        poly[i] = p[i+10];
        poly[i+kk] = p[i+kk+10];
    }
    for (i=0; i<n; i++) {
        yab(x, &i, &kk, poly, &a, &b);  // refine limits here 
        p[6] = x[i];                    // pass forward the value of x; consider &ex etc. 
        Rdqags(fy, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
        x[i] = result;
    }
}

double integral2Dcpp  (const int fn, 
                       const int m, 
                       const int c, 
                       const NumericMatrix &gsbval, 
                       const NumericMatrix &traps,
                       const NumericMatrix &mask, 
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

//===============================================================
// void integral2Dtest
//     (int *fn, int *m, int *c, double *gsbval, int *cc, double *traps,
//     double *mask, int *n1, int *n2, int *kk, int *mm, double *result)
// {
//     double *ex;
//     double ax=1e20;
//     double bx=-1e20;
//     double res;
//     double epsabs = 0.0001;
//     double epsrel = 0.0001;
//     double abserr = 0;
//     int neval = 0;
//     int ier = 0;
//     int limit = 100;
//     int lenw = 400;
//     int last = 0;
//     int iwork[100];
//     double work[400];
//     int k;
//     int ns;
// 
//     // limits from bounding box of this polygon 
//     ns = *n2 - *n1 + 1;
//     for (k=0; k<ns; k++) {
//         ax = R::fmin2(ax, traps[k+ns]);
//         bx = R::fmax2(bx, traps[k+ns]);
//     }
//     ex = (double *) R_alloc(10 + 2 * *kk, sizeof(double));
//     ex[0] = gsbval[*c];   // 1.0? 
//     ex[1] = gsbval[*cc + *c];
//     ex[2] = gsbval[2* *cc + *c];
//     ex[3] = *fn;
//     ex[4] = mask[*m];
//     ex[5] = mask[*m+ *mm];
//     ex[6] = 0;
//     ex[7] = 0;
//     ex[8] = 0;
//     ex[9] = ns;
// 
//     for (k=0; k<ns; k++) {
//         ex[k+10] = traps[k+ *n1];
//         ex[k+ns+10] = traps[k+ *n1 + *kk];
//     }
//     Rdqags(fx, ex, &ax, &bx, &epsabs, &epsrel, &res, &abserr, &neval, &ier,
//           &limit, &lenw, &last, iwork, work);
//     *result = res;
// }

//===============================================================

void fx1 (double *x, int n, void *ex) {
    int i;
    int ns;
    int fn;
    rpoint line[maxvertices * 2];
    rpoint mxy;
    rpoint xy;
    double * p;
    double * cumd;
    double d;
    fnptr fnzr; // = zhnr;
    // extract parameters passed in void pointer ex 
    NumericVector gsb(3);   // passed to fnzr
    p = (double*) ex;
    gsb[0] = p[0];
    gsb[1] = p[1];
    gsb[2] = p[2];
    fn = round(p[3]);
    mxy.x = p[4];
    mxy.y = p[5];
    ns = round(p[9]);
    // coordinates of vertices 
    for (i=0; i<ns; i++) {
        line[i].x = p[i+10];
        line[i].y = p[i+ns+10];
    }
    // cumulative distance along line 
    cumd = (double *) R_alloc(ns + 1, sizeof(double)); 
    cumd[0] = 0;
    for (i=0; i<(ns-1); i++) {
        cumd[i+1] = cumd[i] + distance1 (line[i],line[i+1]);
    }
    // set detection function - default zhnr 
    fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
    // for each x in x[] 
    for (i=0; i<n; i++) {
      xy = getxy (x[i], cumd, line, ns, 0);
      d = distance1 (xy, mxy);
      x[i] = fnzr(gsb, d);   // z(r) 
    }
}

double integral1Dcpp
    (const int fn, 
     const int m, 
     const int c, 
     const NumericMatrix &gsbval, 
     const NumericMatrix &traps,
     const NumericMatrix &mask, 
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

// [[Rcpp::export]]
bool ontransectcpp (
    NumericVector xy,
    NumericMatrix transect,
    int    n1,
    int    n2,
    double tol)
{
  // Is point xy on transect?
  // We assume transect coordinates are in col-major order (x's then y's)
  // http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/ 2009-11-29
  // http://paulbourke.net/geometry/pointlineplane/ 2019-10-13
  int k;
  double r;
  double u;
  rpoint p,p1,p2,p3;
  double minr = 1e20;
  if (n2>=transect.nrow()) stop ("invalid input ontransectcpp");
  
  p3.x = xy(0);
  p3.y = xy(1);
  
  for (k= n1; k < n2; k++)
  {
    p1.x = transect(k,0);
    p1.y = transect(k,1);
    p2.x = transect(k+1,0);
    p2.y = transect(k+1,1);
    if (distance1(p1,p2) > 0) {
      u = ((p3.x-p1.x) * (p2.x-p1.x) + (p3.y-p1.y) * (p2.y-p1.y)) /
        ((p2.x-p1.x) * (p2.x-p1.x) + (p2.y-p1.y) * (p2.y-p1.y));
      if ((u>=0) && (u<=1)) {
        p.x = p1.x + u * (p2.x-p1.x);
        p.y = p1.y + u * (p2.y-p1.y);
        r = distance1 (p,p3);
        minr = R::fmin2(r,minr);
      }
    }
  }
  // check at each vertex to be sure 
  for (k= n1; k <= n2; k++) {
    p1.x = transect(k,0);
    p1.y = transect(k,1);
    r = distance1 (p1,p3);
    minr = R::fmin2(r,minr);
  }
  return(minr < tol);
}
//----------------------------------------------------------------

// [[Rcpp::export]]
double alongtransectcpp (
    NumericVector xy,
    NumericMatrix transect,
    int    n1,
    int    n2,
    double tol
    )
{
  // How far is point xy from start of transect? 2011-06-07
  // We assume transect coordinates are in col-major order (x's then y's)
  // http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/ 29/11/09
  int k;
  double r;
  double u;
  rpoint p,p1,p2,p3;
  
  double along=0.0;
  if (n2>=transect.nrow()) stop ("invalid input alongtransectcpp");
  p3.x = xy(0);
  p3.y = xy(1);
  
  for (k= n1; k < n2; k++)
  {
    p1.x = transect(k,0);
    p1.y = transect(k,1);
    r = distance1 (p1,p3);
    if (r < tol) {
      return(along);
    } 
    p2.x = transect(k+1,0);
    p2.y = transect(k+1,1);
    if (distance1(p1,p2) > 0) {
      u = ((p3.x-p1.x) * (p2.x-p1.x) + (p3.y-p1.y) * (p2.y-p1.y)) /
        ((p2.x-p1.x) * (p2.x-p1.x) + (p2.y-p1.y) * (p2.y-p1.y));
      if ((u>=0) && (u<=1)) {
        p.x = p1.x + u * (p2.x-p1.x);
        p.y = p1.y + u * (p2.y-p1.y);
        r = distance1 (p,p3);
        if (r < tol) {
          along += distance1(p,p1);
          return(along);
        } 
      }
      along += distance1(p1,p2);
    }
  }
  return(along);
}
//----------------------------------------------------------------

void rgr(double *x, int n, void *ex) {
  int i;
  int fn;
  double * p;
  NumericVector tmp(4);
  p = (double*) ex;
  for (i=0; i<4; i++) tmp[i] = p[i];
  fn = tmp[3];
  fnptr fnp; // = zhnr;
  fnp = getzfnr(fn);
  for (i=0; i<n; i++) {
    x[i] = x[i] * fnp(tmp,x[i]);   // r.h(r) 
  }
}

//===============================================================
double hintegralcpp (
    const int fn, 
    const NumericVector &gsb) 
  // integral of radial 2-D function 
  // from 2017-03-22 this is strictly a hazard function 
{
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
  b = 1;    // signals bounds 0,Inf 
  ex[0] = gsb(0);
  ex[1] = gsb(1);
  ex[2] = gsb(2);
  ex[3] = fn;
  
  Rdqagi(rgr, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
         &limit, &lenw, &last, iwork, work);
  
  // ignoring ier etc. 
  return (result * 2 * M_PI);
}
//----------------------------------------------------------------


void justgr(double *x, int n, void *ex) {
  int i;
  int fn;
  double * p;
  NumericVector tmp(4);
  p = (double*) ex;
  for (i=0; i<4; i++) tmp[i] = p[i];
  fn = tmp[3];
  fnptr fnp; // = zhnr;
  fnp = getzfnr(fn);
  for (i=0; i<n; i++) {
    x[i] = fnp(tmp,x[i]);   // h(r) 
  }
}

// integral of 1-D function 
// from 2017-03-22 this is strictly a hazard function 
double hintegral1cpp (
    const int fn, 
    const NumericVector &gsb) 
{
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
  
  // treat uniform separately 
  if (fn == 4) {
    return (2 * gsb(0) * gsb(1));  // symmetric about 0 
  }
  
  a = 0;
  b = 1;    // signals bounds 0,Inf 
  ex[0] = gsb(0);
  ex[1] = gsb(1);
  ex[2] = gsb(2);
  ex[3] = fn;
  
  Rdqagi(justgr, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, 
         &ier, &limit, &lenw, &last, iwork, work);
  
  // ignoring ier etc. 
  return (result * 2);    // symmetric about 0 
}
//----------------------------------------------------------------
