// include guard
#ifndef __secr_h_INCLUDED__   // if secr.h hasn't been included yet...
#define __secr_h_INCLUDED__   // #define this so the compiler knows it has been included

//------------------------------------------------------------------------------
// BOOST used for statistical distributions from secr 4.4.6 October 2021
// return NAN for invalid inputs
// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html
#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
// must follow define domain error policy...
#include <boost/math/distributions.hpp>     
//------------------------------------------------------------------------------

#include <Rcpp.h>
#include <RcppParallel.h>

// #include <R.h>       // random numbers 

using namespace Rcpp;
using namespace RcppParallel;

// constants
#define fuzz 1e-200
#define huge 1e10
#define maxnpoly 1000   
#define maxnmix 2    
#define maxvertices 400

//-------------------
// data structures   
//-------------------
struct trap_animal {
    int     trap;
    int     animal;
    double  time;
};
struct rpoint {
    double x;
    double y;
};

//--------------------------------------------------------------------------

int i3 (int i, int j, int k, int ii, int jj);
int i4 (int i, int j, int k, int l, int ii, int jj, int kk);

//------------------------------------------------------
// detectfn.cpp 
//------------------------------------------------------
typedef double (*fnptr)(const Rcpp::NumericVector&, const double);
typedef double (*fnptrC)(const std::vector<double>&, const double);
fnptr getzfnr (int fn);
fnptrC getgfns (int fn);
fnptrC getzfnrC (int fn);

double pfnS (
        const int fn,
        const double d2val,
        const std::vector<double> &gsb,
        const std::vector<double> &miscparm,
        const double w2);
    
//--------------------------------------------------------------------------

// not used 2021-10-17
// customised dnbinom parameterised as size, mu 
// double gnbinom (int count, int size, double mu, int uselog);
//--------------------------------------------------------------------------

// probability of count with distribution specified by binomN 
double countp (int count, int binomN, double lambda);
//--------------------------------------------------------------------------

double mufn (
    const int k,
    const int m,
    const double b0,
    const double b1,
    const Rcpp::NumericMatrix &A1,
    const Rcpp::NumericMatrix &A2,
    const int spherical);

//---------------------------------------------------------------------

double mufnL (
    const int k,
    const int m,
    const double b0,
    const double b1,
    const Rcpp::NumericMatrix &dist2,
    const int spherical);

//---------------------------------------------------------------------

double SegCircle2 (double p1x, double p1y, double p2x, double p2y, double scx, double scy, 
    double r);

//---------------------------------------------------------------------

double expmin (double x);

double distance1 (const rpoint p1, const rpoint p2);

double rcount (const int binomN, const double lambda, const double Tsk);
//---------------------------------------------------------------------

rpoint getxy(
    const double l, 
    double cumd[], 
    const rpoint line[], 
    const int kk, 
    const int offset);   // double before 2022-01-18

//---------------------------------------------------------------------

double randomtime (double p);
double randomtimel (double lambda);
//---------------------------------------------------------------------

void probsort (
    const int n, 
    std::vector<trap_animal> &tran);

//---------------------------------------------------------------------

double gr (
    const int fn, 
    Rcpp::NumericVector gsb, 
    const rpoint xy, 
    const rpoint animal);
//---------------------------------------------------------------------

// random point from 2-D radial distribution specified by g function 
Rcpp::NumericVector gxy (const int fn, 
    const Rcpp::NumericVector par, 
    const double w);
  
//---------------------------------------------------------------------

double hazard (double pp);

//---------------------------------------------------------------------
 
double gpois (int count, double lambda);
double gbinom(int count, int size, double p);
double pski ( int binomN, int count, double Tski, double g, double pI);

//--------------------------------------------------------------------------

// double d2 (int k, int m, double A1[], double A2[], int A1rows, int A2rows);

double d2cpp (
    const int k, 
    const int m, 
    const Rcpp::NumericMatrix &A1, 
    const Rcpp::NumericMatrix &A2);

Rcpp::List makelookupcpp (
    const Rcpp::NumericMatrix &x);  


// Functions to characterize detector type 
// polygon, transect and signal detector types must be constant across occasions

bool anypolygon (const Rcpp::IntegerVector detect);
bool anytransect (const Rcpp::IntegerVector detect);
// bool anysignal (const Rcpp::IntegerVector detect);
bool anytelemetry (const Rcpp::IntegerVector detect);
bool allpoint (const Rcpp::IntegerVector detect, bool allowsignal, bool allowtelem);

// miscellaneous functions

bool insidecpp (
        const Rcpp::NumericVector &xy,
        const int    n1,
        const int    n2,
        const Rcpp::NumericMatrix &poly);

void fillngcpp(const int nc, 
               const int gg, 
               const Rcpp::IntegerVector &grp, 
               std::vector<int> &ng);

void yab(double x[], int *i, int *np, double poly[], double *a, double *b);
void fy(double *x, int n, void *ex);
void fx(double *x, int n, void *ex);
void fx1 (double *x, int n, void *ex);
    
#endif  // __secr_h_INCLUDED__
