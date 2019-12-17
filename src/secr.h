// include guard
#ifndef __secr_h_INCLUDED__   // if secr.h hasn't been included yet...
#define __secr_h_INCLUDED__   // #define this so the compiler knows it has been included

#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;

#include <R.h>       // random numbers 
#include <Rmath.h>   // R math functions e.g. dbinom, dpois 

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
typedef double (*fnptr)(const Rcpp::NumericVector, const double);
typedef double (*fnptrC)(const std::vector<double>, const double);
fnptr getgfnr (int fn);
fnptr getzfnr (int fn);
fnptrC getzfnrC (int fn);
double pfn (
        const int fn, 
        const double d2val, 
        const Rcpp::NumericVector &gsb,
        const Rcpp::NumericVector &miscparm, 
        const double w2);

int par3 (int fn);

//--------------------------------------------------------------------------

// customised dnbinom parameterised as size, mu 
double gnbinom (int count, int size, double mu, int uselog);
//--------------------------------------------------------------------------

// binomial density allowing non-integer (floating point) size 
// double gbinomFP (int count, double size, double p, int uselog);
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

void SegCircle (double *p1x, double *p1y, double *p2x, double *p2y, double *scx, double *scy, 
    double *r, double *seg);
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
    const double offset);

rpoint getxycpp(
        const double l, 
        const std::vector<double> &cumd, 
        const RcppParallel::RMatrix<double> &line, 
        const int kk, 
        const double offset);
//---------------------------------------------------------------------

double randomtime (double p);
//---------------------------------------------------------------------

void probsort (
    const int n, 
    std::vector<trap_animal> &tran);

//---------------------------------------------------------------------

double gr (
    const int fn, Rcpp::NumericVector, 
    const rpoint xy, 
    const rpoint animal);
//---------------------------------------------------------------------

// random point from 2-D radial distribution specified by g function 
NumericVector gxy (const int fn, 
    const Rcpp::NumericVector par, 
    const double w);
  
//---------------------------------------------------------------------

double hazard (double pp);

void getdetspec (
    const Rcpp::IntegerVector &detect, 
    const int fn, 
    const int nc,  
    const int nc1, 
    const int cc, 
    const int nmix, 
    const int nd, 
    const int nk, 
    const int ss, 
    const int mm, 
    const Rcpp::IntegerVector &PIA, 
    const Rcpp::NumericVector &miscparm, 
    const std::vector<int> &start, 
    std::vector<double> &detspec);

void geth2 (
    const int nc1, 
    const int cc, 
    const int nmix, 
    const int mm, 
    const Rcpp::IntegerVector &PIA, 
    const std::vector<double> &hk, 
    const Rcpp::NumericMatrix &Tsk, 
    std::vector<double> &h, 
    std::vector<int> &hindex);

// int getstart(
//     const Rcpp::IntegerVector &detect, 
//     std::vector<int> &start, 
//     const int nc1, 
//     const int nc, 
//     const int ss, 
//     const int nk, 
//     const Rcpp::IntegerVector &w);

//---------------------------------------------------------------------
// 
double gpois (int count, double lambda, int uselog);
double gbinom(int count, int size, double p, int uselog);
double gbinomFP (int count, double size, double p, int uselog);
double pski ( int binomN, int count, double Tski, double g, double pI);

//--------------------------------------------------------------------------

// double d2 (int k, int m, double A1[], double A2[], int A1rows, int A2rows);

double d2cpp (
    const int k, 
    const int m, 
    const NumericMatrix &A1, 
    const NumericMatrix &A2);

List makelookupcpp (
    const NumericMatrix &x);  


// Functions to characterize detector type 
// polygon, transect and signal detector types must be constant across occasions

bool anyexclusive (const IntegerVector detect);
bool anycapped (const IntegerVector detect);
bool anypolygon (const IntegerVector detect);
bool anytransect (const IntegerVector detect);
bool anysignal (const IntegerVector detect);
bool anytelemetry (const IntegerVector detect);
bool alltelemetry (const IntegerVector detect);
bool allpobool (const IntegerVector detect, bool allowsignal, bool allowtelem);
bool allcapped  (const IntegerVector detect);
bool allmulti (const IntegerVector detect);
bool allpoint (const IntegerVector detect, bool allowsignal, bool allowtelem);

bool anyvarying (const int nc, const int ss, const int nk, const int nmix,
                 const IntegerVector &PIA0);
bool anyb (const NumericMatrix &gsbval, const NumericMatrix &gsb0val);

// miscellaneous functions

// std::vector<int> fillcumkcpp(
//     const IntegerVector detect, 
//     const int ss, 
//     const IntegerVector kk);
  
int nval(int detect0, int nc1, int cc, int ss, int nk);

NumericMatrix makedist2cpp (
    const NumericMatrix &traps, 
    const NumericMatrix &mask);

void squaredistcpp (NumericMatrix &dist2);

bool insidecpp (
        const NumericVector &xy,
        const int    n1,
        const int    n2,
        const Rcpp::NumericMatrix &poly);

void fillngcpp(const int nc, 
               const int gg, 
               const IntegerVector &grp, 
               std::vector<int> &ng);

//---------------------------------------------------------------
// Return probability individual n belongs to class x. This may be binary 
//   (0/1) in the case of known class, or continuous if class is unknown 
double classmembership (
    const int n, 
    const int x, 
    const IntegerVector &knownclass, 
    const std::vector<double> &pmixn, 
    const int nmix);
  
void yab(double x[], int *i, int *np, double poly[], double *a, double *b);
void fy(double *x, int n, void *ex);
void fx(double *x, int n, void *ex);
void fx1 (double *x, int n, void *ex);
    
#endif  // __secr_h_INCLUDED__
