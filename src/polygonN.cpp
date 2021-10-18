#include "poly.h"
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;
using namespace Numer;

//#include <R_ext/Utils.h>
//#include <R_ext/Applic.h>


// C++ code for R package 'secr' 
// Murray Efford 

// This file contains functions for integration of home range overlap
// with polygon detectors and similar.

// Also exports utility functions ontransectcpp and alongtransectcpp

//---------------------------------------------------------------------------------
// original 2-D code using Rdqags
//---------------------------------------------------------------------------------

// find upper and lower points on perimeter of poly at x-coordinate x
// void yab(double x[], int *i, int *np, double poly[], double *a, double *b) {
//     int k;
//     int nv = 0;
//     double ab[3];
//     // note 'sign' is RMath function
//     for (k=0; k< (*np-1); k++) {
//         if (R::sign(poly[k]- x[*i]) != R::sign(poly[k+1]- x[*i])) {
//             ab[nv] = poly[k+ *np] + (x[*i]-poly[k]) * (poly[k+1+*np]-poly[k+*np]) /
//                 (poly[k+1]-poly[k]);
//             nv++;
//         }
//         if (nv>2) break;
//     }
//     if (ab[0]>ab[1])
//     {*a = ab[1]; *b = ab[0];}
//     else
//     {*a = ab[0]; *b = ab[1];}
// 
// }
// 
// void fy(double *x, int n, void *ex) {
//     int i;
//     int fn;
//     double * p;
//     double mx,my;
//     double xy[2];
//     double d;
//     fnptr fnzr; // = zhnr;
//     NumericVector gsb(3);   // passed to fnzr
//     p = (double*) ex;
//     gsb[0] = p[0];
//     gsb[1] = p[1];
//     gsb[2] = p[2];
//     fn = round(p[3]);
//     mx = p[4];
//     my = p[5];
//     xy[0] = p[6];
//     // set detection function
//     fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
//     for (i=0; i<n; i++) {
//         xy[1] = x[i];   // set each y value
//         d = std::sqrt ( (xy[1]-my)*(xy[1]-my) + (xy[0]-mx)*(xy[0]-mx) );
//         x[i] = fnzr(gsb, d);   // z(r)
//     }
// }
// 
// void fx(double *x, int n, void *ex) {
//     int i;
//     double * p;
// 
//     double poly[maxvertices * 2];
//     double a;
//     double b;
// 
//     double epsabs = 0.0001;
//     double epsrel = 0.0001;
//     double result = 0;
//     double abserr = 0;
//     int neval = 0;
//     int ier = 0;
//     int limit = 100;
//     int lenw = 400;
//     int last = 0;
//     int iwork[100];
//     double work[400];
//     int kk;
//     p = (double*) ex;
//     kk = round(p[9]);
// 
//     for (i=0; i<kk; i++) {
//         poly[i] = p[i+10];
//         poly[i+kk] = p[i+kk+10];
//     }
//     for (i=0; i<n; i++) {
//         yab(x, &i, &kk, poly, &a, &b);  // refine limits here
//         p[6] = x[i];                    // pass forward the value of x; consider &ex etc.
//         Rdqags(fy, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
//                &limit, &lenw, &last, iwork, work);
//         x[i] = result;
//     }
// }
// 
// double integral2Dcpp  (const int &fn,
//                        const int &m,
//                        const int &c,
//                        const RMatrix<double> &gsbval,
//                        const RMatrix<double> &traps,
//                        const RMatrix<double> &mask,
//                        const int &n1,
//                        const int &n2,
//                        double ex[]) {
//     double ax=1e20;
//     double bx=-1e20;
//     double epsabs = 0.0001;
//     double epsrel = 0.0001;
//     double result = 0;
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
//     int reportier = 0;
// 
//     // limits from bounding box of this polygon
//     ns = n2-n1+1;
//     for (k=0; k<ns; k++) {
//         ax = R::fmin2(ax, traps[k+n1]);
//         bx = R::fmax2(bx, traps[k+n1]);
//     }
// 
//     // pass parameters etc. through pointer
//     ex[0] = gsbval(c,0);
//     ex[1] = gsbval(c,1);
//     ex[2] = gsbval(c,2);
//     ex[3] = fn;
//     ex[4] = mask(m,0);
//     ex[5] = mask(m,1);
//     ex[6] = 0;
//     ex[7] = 0;
//     ex[8] = 0;
//     ex[9] = ns;
// 
//     // also pass polygon vertices
//     for (k=0; k<ns; k++) {
//         ex[k+10] = traps(k+n1,0);        // x
//         ex[k+ns+10] = traps(k+n1,1);       // y
//     }
//     Rdqags(fx, ex, &ax, &bx, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
//            &limit, &lenw, &last, iwork, work);
//     if ((ier != 0) & (reportier))
//         Rprintf("ier error code in integral2D %5d\n", ier);
//     return (result);
// }
// // end original 2-D code using Rdqags
// //---------------------------------------------------------------------------------

bool insidecppC (
    const Numer::Constvec &xy,
    const int    &n1,
    const int    &n2,
    const RcppParallel::RMatrix<double> &poly)
{
  // Is point xy inside poly?
  // Based on contribution on s-news list by Peter Perkins 23/7/96
  // We assume poly is closed, and in col-major order (x's then y's)
  
  double theta = 0;
  double cutoff = 1e-6;
  int k;
  int ns;
  double N;
  double d;
  ns = n2 - n1 + 1;   // number of selected points 
  std::vector<double> temp((ns+1) * 2);
  
  // get & translate to coords centered at each test point 
  for (k=0; k < ns; k++)
  {
    temp[k]      = poly(k + n1,0) - xy[0];    // x 
    temp[k + ns] = poly(k + n1,1) - xy[1];    // y 
  }
  
  for (k=0; k < (ns-1); k++)
  {
    N = temp[k] * temp[k+1 + ns] - temp[k + ns] * temp[k+1];
    d = temp[k] * temp[k+1]      + temp[k + ns] * temp[k+1 + ns];
    if (fabs(d)>0) { N = N/fabs(d);  d = d/fabs(d); }
    theta += std::atan2(N, d);
  }
  theta = fabs(theta);
  return (fabs(theta - 2* M_PI) < cutoff);    // M_PI is cmath.h constant 
}
//--------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// 2-D code using repeated 1-D RcppNumerical Numer::integrate Func
//---------------------------------------------------------------------------------

class yslice: public Func
{
private:
    std::vector<double> gsb;
    int fn;
    double mx, my;
    double x;
    fnptrC fnzr; // = zhnr;

public:
    yslice(const std::vector<double> gsb,
           const int fn,
           const double mx,
           const double my,
           const double x)
        : gsb(gsb), fn(fn), mx(mx), my(my), x(x) {
        // set detection function
        fnzr = getzfnrC(fn);
    }

    double operator()(const double& y) const
    {
        double d;
        d = std::sqrt ((y-my)*(y-my) + (x-mx)*(x-mx));
        return(fnzr(gsb, d));   // z(r)
    }
};

class xfn: public Func {
    
private:
    std::vector<double> gsb;
    RMatrix<double> poly;
    const int n1;
    const int n2;
    int fn;
    double mx;
    double my;
    fnptrC fnzr; // = zhnr;
    int np;
    
    class yslicei: public Func
    {
    private:
        const std::vector<double> gsb;
        const int fn;
        int n1;
        int n2;
        const double mx;
        const double my;
        fnptrC fnzr; // = zhnr;
        
    public:
        
        double x=0.0;
        
        yslicei(const std::vector<double> &gsb,
                const int &fn,
                const double &mx,
                const double &my): gsb(gsb), fn(fn), mx(mx), my(my)
        {fnzr = getzfnrC(fn);}
        
        double operator()(const double& y) const
        {
            double d;
            d = std::sqrt (pow(y-my,2) + pow(x-mx,2));
            return(fnzr(gsb, d));   // z(r)
        }
    };
    
public:
    
    xfn(const std::vector<double> &gsb,
        const RMatrix<double> &poly,
        const int &n1,
        const int &n2,
        const int &fn,
        const double &mx,
        const double &my
    )
        : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my) {
        // set detection function
        fnzr = getzfnrC(fn);
        np = poly.nrow();
    }
    
    
    // find upper and lower points on perimeter of poly at x-coordinate x
    std::vector<double> ylim(const double &x) const {
        int nv = 0;
        std::vector<double> ab(2);
        double tmp;
        for (int k=n1; k<n2; k++) {
            // note 'sign' is R function best avoided
            // if (R::sign(poly[k]- x) != R::sign(poly[k+1]- x)) {
            if (((poly[k]< x) && (poly[k+1] > x)) || ((poly[k]> x) && (poly[k+1] < x))) {
                ab[nv] = poly[k+ np] + (x-poly[k]) * (poly[k+1+np]-poly[k+np]) /
                    (poly[k+1]-poly[k]);
                nv++;
            }
            if (nv>2) break;
        }
        if (ab[0]>ab[1])
        {
            tmp = ab[1];
            ab[1] = ab[0];
            ab[0] = tmp;
        }
        return(ab);
    }
    
    double operator()(const double &x) const
    {
        double err_est;
        int err_code;
        std::vector<double> lim;
        yslicei f(gsb, fn, mx, my);
        f.x = x;
        lim = ylim(x);  // refine limits here
        const double res = integrate(f, lim[0], lim[1], err_est, err_code);
        return res;
    }
};

//---------------------------------------------------------------------------------------
// xfn2 is a version of xfn that uses pointwise insidecppC rather than merely integrating
// between upper and lower bounds in y dimension

class xfn2: public Func {
    
private:
    std::vector<double> gsb;
    RMatrix<double> poly;
    int n1;
    int n2;
    int fn;
    double mx;
    double my;
    double ay;
    double by;
    
    fnptrC fnzr; // = zhnr;
    int np;
    
    class yslicei: public Func
    {
    private:
        const std::vector<double> gsb;
        const RMatrix<double> &poly;
        const int n1;
        const int n2;
        const int fn;
        const double mx;
        const double my;
        const double ay;
        const double by;
        fnptrC fnzr; // = zhnr;
        
    public:
        
        double x=0.0;
        
        yslicei(const std::vector<double> &gsb,
                const RMatrix<double> &poly,
                const int &n1,
                const int &n2,
                const int &fn,
                const double &mx,
                const double &my,
                const double &ay,
                const double &by)
            : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my), ay(ay), by(by)
        {fnzr = getzfnrC(fn);}
        
        double operator()(const double& y) const
        {
            double d;
            Eigen::VectorXd xy(2);
            xy << x,y;
            if (insidecppC(xy, n1, n2, poly)) {
                d = std::sqrt (pow(y-my,2) + pow(x-mx,2));
                return(fnzr(gsb, d));   // z(r)
            }
            else {
                return(0);
            }
        }
    };
    
public:
    
    xfn2(const std::vector<double> &gsb,
        const RMatrix<double> &poly,
        const int &n1,
        const int &n2,
        const int &fn,
        const double &mx,
        const double &my,
        const double &ay,
        const double &by
    )
        : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my), ay(ay), by(by) {
        // set detection function
        fnzr = getzfnrC(fn);
        np = poly.nrow();
    }

    double operator()(const double &x) const
    {
        double err_est;
        int err_code;
        std::vector<double> lim;
        yslicei f(gsb, poly, n1, n2, fn, mx, my, ay, by);
        f.x = x;
        const double res = integrate(f, ay, by, err_est, err_code);
        return res;
    }
};
//---------------------------------------------------------------------------------------

// Apply either xfn or xfn2

double integral2DNRcpp
    (const int &fn,
     const int &m,
     const int &c,
     const RMatrix<double> &gsbval,
     const RMatrix<double> &poly,
     const RMatrix<double> &mask,
     const int &n1,
     const int &n2,
     const bool &convex)
{
    double res;
    double err_est;
    int err_code;
    double ax=1e10;
    double bx=-1e10;
    double ay=1e10;    // used only with xfn2
    double by=-1e10;   // used only with xfn2
    int ns = n2-n1+1;
    for (int i=0; i<ns; i++) {
        ax = std::min(ax, poly(n1+i,0));
        bx = std::max(bx, poly(n1+i,0));
        ay = std::min(ay, poly(n1+i,1));    // used only with xfn2
        by = std::max(by, poly(n1+i,1));    // used only with xfn2
    }

    std::vector<double> gsb(4);
    int npar = gsbval.ncol();
    for (int i=0; i<npar; i++) gsb[i] = gsbval(c,i);

    // "Possible values are GaussKronrod{15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201}"
    // GaussKronrod41 default
    // GaussKronrod15 causes R failure
    // this variation seems slower
    // res = integrate(f, ax, bx, err_est, err_code, 100, 1e-6, 1e-5, Integrator<double>::GaussKronrod21);
    if (convex) {
        xfn f(gsb, poly, n1, n2, fn, mask(m,0), mask(m,1));
        res = integrate(f, ax, bx, err_est, err_code);
    }
    else {
        xfn2 f2(gsb, poly, n1, n2, fn, mask(m,0), mask(m,1), ay, by);
        res = integrate(f2, ax, bx, err_est, err_code);
    }
    return (res);
}

//===============================================================

class fx1func: public Numer::Func
{
private:
    std::vector<double> gsb;      
    RMatrix<double> line;
    int fn;
    double mx, my;
    fnptrC fnzr;
    int ns;
    std::vector<double> cumd; 
public:
    std::vector<double> cumdistance(RMatrix<double> line) {
        int ns = line.nrow();
        std::vector<double> cumd (ns,0);
        for (int i=0; i<(ns-1); i++) {
            cumd[i+1] = cumd[i] + 
                std::sqrt (pow(line(i,0)-line(i+1,0),2) + 
                pow(line(i,1)-line(i+1,1),2));
        }
        return(cumd);
    }
    
    fx1func(const std::vector<double> &gsb,
            const RMatrix<double> &line,
            const int fn,
            const double mx,
            const double my
    ) 
        : gsb(gsb), line(line), fn(fn), mx(mx), my(my) {
        fnzr = getzfnrC(fn);   // set detection function 
        ns = line.nrow();      // cumulative distance along line 
        cumd = cumdistance(line);
    }
    
    
    double operator()(const double& x) const
    {
        double d;
        rpoint xy;
        xy = getxycpp (x, cumd, line, ns, 0);
        d = std::sqrt (pow(xy.x-mx,2) + pow(xy.y-my,2));
        return(fnzr(gsb, d));   // z(r) 
    }
};

double integral1DNRcpp
    (const int fn, 
     const int m, 
     const int c, 
     const RMatrix<double> &gsbval, 
     const RMatrix<double> &traps,
     const RMatrix<double> &mask, 
     const int n1, 
     const int n2)
{
    double err_est;
    int err_code;
    double ax=0;
    double bx=0;
    int k;
    
    if (gsbval.ncol()>4) Rcpp::stop("bad gsbval matrix");
    std::vector<double> gsb(4);
    // uniform treated separately 
    if (fn == 4) {
        for (k=n1+1; k<=n2; k++) {  // upper bound is length of this transect 
            bx += SegCircle2(traps(k-1,0), traps(k-1,1), traps(k,0), traps(k,1),
                             mask(m,0), mask(m,1), gsbval(c,1));
        }
        return (bx);
    }
    
    for (k=n1+1; k<=n2; k++) {  // upper bound is length of this transect 
        bx += std::sqrt( (traps(k,0) - traps(k-1,0)) * (traps(k,0) - traps(k-1,0)) +
            (traps(k,1) - traps(k-1,1)) * (traps(k,1) - traps(k-1,1)) );
    }
    int npar = gsbval.ncol();
    for (int i=0; i<npar; i++) gsb[i] = gsbval(c,i);
    fx1func f(gsb, traps, fn, mask(m,0), mask(m,1));
    const double res = integrate(f, ax, bx, err_est, err_code);
    return (res);
}


//---------------------------------------------------------------------------------
// alternative using 2-D code using RcppNumerical Numer::integrate MFunc
//---------------------------------------------------------------------------------

class fx2func: public Numer::MFunc
{
private:
    std::vector<double> gsb;
    RMatrix<double> poly;
    int n1;
    int n2;
    int fn;
    double mx;
    double my;
    fnptrC fnzr;
    
public:
    fx2func(const std::vector<double> &gsb,
            const RMatrix<double> &poly,
            const int &n1,
            const int &n2,
            const int &fn,
            const double &mx,
            const double &my
    ) 
        : gsb(gsb), poly(poly), n1(n1), n2(n2), fn(fn), mx(mx), my(my) {
        fnzr = getzfnrC(fn);    // set detection function 
    }
    
    double operator()(Numer::Constvec &x) const
    {
        double d;
        if (insidecppC(x, n1, n2, poly)) {
            d = std::sqrt (pow(x[0]-mx,2) + pow(x[1]-my,2));
            return(fnzr(gsb, d));   // z(r) 
        }
        else
            return 0;
    }
};

// double integral2DNRcpp  (
//         const int &fn,
//         const int &m,
//         const int &c,
//         const RMatrix<double> &gsbval,
//         const RMatrix<double> &poly,
//         const RMatrix<double> &mask,
//         const int &n1,
//         const int &n2) 
// {
//     Eigen::VectorXd lower(2);
//     Eigen::VectorXd upper(2);
//     double err_est;
//     int err_code;
//     int k;
//     int ns;
//     double xmin,ymin = 1e100;
//     double xmax,ymax = -1e100;
// 
//     // limits from bounding box of this polygon
//     ns = n2-n1+1;
//     for (k=0; k<ns; k++) {
//         xmin = std::min(xmin, poly(k+n1,0));  
//         ymin = std::min(ymin, poly(k+n1,1));  
//         xmax = std::max(xmax, poly(k+n1,0));  
//         ymax = std::max(ymax, poly(k+n1,1));  
//     }
//     lower[0] = xmin;
//     lower[1] = ymin;
//     upper[0] = xmax;
//     upper[1] = ymax;
//     
//     std::vector<double> gsb(4,0);
//     int npar = gsbval.ncol();
//     for (int i=0; i<npar; i++) gsb[i] = gsbval(c,i);
//     
//     fx2func f(gsb, poly, n1, n2, fn, mask(m,0), mask(m,1));
//     
//     const double res = integrate(f, lower, upper, err_est, err_code);
//     return (res);
// }
//===============================================================

// Miscellaneous 1-D integrations

class rgrfn: public Func
{
private:
    std::vector<double> gsb;
    int fn;
    fnptrC fnzr; // = zhnr;
    
public:
    rgrfn(const std::vector<double> gsb,
          const int fn) 
        : gsb(gsb), fn(fn) {
        // set detection function 
        fnzr = getzfnrC(fn);  
    }
    double operator()(const double& r) const
    {
        return( r * fnzr(gsb, r));   // z(r) 
    }
};

double hintegral2Ncpp (
    const int fn, 
    const std::vector<double> &gsb) 
  // integral of radial 2-D function 
{
  double err_est;
  int err_code;
  rgrfn f(gsb, fn);
  const double res = integrate(f, 0, 20*gsb[1], err_est, err_code);
  return (res * 2 * M_PI);
}
//===============================================================

class grfn: public Func
{
private:
    std::vector<double> gsb;
    int fn;
    fnptrC fnzr; // = zhnr;
    
public:
    grfn(const std::vector<double> gsb,
         const int fn) 
    : gsb(gsb), fn(fn) {
        fnzr = getzfnrC(fn);  // set detection function 
    }
    double operator()(const double& x) const
    {
        return fnzr(gsb, x);   // z(r) 
    }
};

// integral of 1-D function (radial 2-D function)
double hintegral1Ncpp (
        const int fn, 
        const std::vector<double> &gsb) 
{
    double err_est;
    int err_code;
    grfn f(gsb, fn);
    // infinite limit fails
    // but see https://www.r-bloggers.com/numerical-integration-over-an-infinite-interval-in-rcpp-2/
    // const double res = integrate(f, 0, R_PosInf, err_est, err_code);
    // ASSUME gsb[1] is scale
    const double res = integrate(f, 0, 20 * gsb[1], err_est, err_code);
    // Rprintf("res %10.7g err_est %10.7g err_code %4d \n", res, err_est, err_code);
    if (err_code>0) Rcpp::stop ("err_code>0 in hintegral1Ncpp");
    return (res * 2);
}
//===============================================================

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
    if (n2>=transect.nrow()) Rcpp::stop ("invalid input ontransectcpp");
    
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
    if (n2>=transect.nrow()) Rcpp::stop ("invalid input alongtransectcpp");
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
