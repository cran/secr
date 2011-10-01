#include "secr.h"

/* unmarked and presence likelihoods */
/* 2011-10-01 */

/* temporary unmarked likelihood */
/* only constant model *cc == 1 */
/* likelihood assumes indepndence between sites AND between successive occasions,
   but latter limitation should be fixable */

void unmarkedloglik (
    int    *w,           /* capture histories (1:nc, 1:ss, 1:kk) */
    int    *nc,          /* number of rows in w */
    int    *ss,          /* number of occasions */
    int    *kk,          /* number of detectors */
    double *D,           /* Parameter value - density */
    double *g0,          /* Parameter value - p */
    double *sigma,       /* Parameter value - radius */
    double *z,           /* Parameter value - shape (hazard rate, cumulative gamma etc) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential etc.*/
    int    *type,        /* code not used */
    double *value,       /* return value */
    int    *resultcode   /* 0 if OK */
    ) {

    int    k,n,s;
    int    *nsk;
    double integral;
    double lambda;
    double par[3];

    /*===============================================================*/

    *resultcode = 1;  /* generic failure code */
    *value = 0;

    if (*ss > 1) 
        warning ("'unmarked' likelihood does not apply for > 1 occasion");

    /* summarise input w as nsk matrix */
    nsk = (int *) R_alloc(*ss * *kk, sizeof(int));
    for (k = 0; k < *kk; k++) {
        for (s = 0; s < *ss; s++) {
            nsk[k* *ss + s] = 0;
            for (n = 0; n < *nc; n++) {
                nsk[k* *ss + s] += w[i3(n,s,k,*nc,*ss)];
	    }
	}
    }
    par[0] = *g0;
    par[1] = *sigma;
    par[2] = *z;
    integral = gintegral(*fn, par);     /* one occasion */
    for (k = 0; k < *kk; k++) {
        for (s = 0; s < *ss; s++) {
            lambda = *D * integral / 10000;
            *value += dpois(nsk[k * *ss + s], lambda, 1);
        }
    }
    *resultcode = 0;   /* successful termination unmarkedloglik */
}
/*==============================================================================*/

/*
Likelihood function for presence detector type, assuming independence
'simple' version is equivalent to Royle & Nichols 2003 Poisson
'integrated' alows for non-step detection function
MGE 2011-09-30
*/

void rgr2(double *x, int n, void *ex) {
  /* r . (1 - (1-g(r))^(*ss-j)) */
    int i;
    int fn;
    double * p;
    double tmp[5];
    p = (double*) ex;
    for (i=0; i<5; i++) tmp[i] = p[i];
    fn = tmp[3];
    fnptr fnp = hn;
    fnp = gethfn(fn);
    for (i=0; i<n; i++) {
        x[i] = x[i] * (1 - pow(1 - fnp(tmp,x[i]), tmp[4])); 
    }
}
double gintegral2 (int fn, int sj, double par[]) {
/* integral of radial 2-D function rgr2 */
    double ex[5];
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
    ex[4] = sj;
    Rdqagi(rgr2, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
    /* ignoring ier etc. */
    return (result * 2 * M_PI);
}

void presenceloglik (
    int    *w,           /* capture histories (1:nc, 1:ss, 1:kk) */
    int    *nc,          /* number of rows in w */
    int    *ss,          /* number of occasions */
    int    *kk,          /* number of detectors */
    double *D,           /* Parameter value - density */
    double *g0,          /* Parameter value - p */
    double *sigma,       /* Parameter value - radius */
    double *z,           /* Parameter value - shape (hazard rate, cumulative gamma etc) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential etc.*/
    int    *type,        /* code 0 = simple, 1 = integrated */
    double *value,       /* return value */
    int    *resultcode   /* 0 if OK */
)

{
    int    k,j,n,s;
    int    *y;           /* vector of detector-specific counts, length *kk 0 <= yk <= ss) */
    double tempsum;
    double lambda;
    double par[3];
    double mu;
    double integral;

    /*===============================================================*/


    /* summarise input w as y vector */
    y = (int *) R_alloc(*kk, sizeof(int));
    for (k = 0; k < *kk; k++) {
        y[k] = 0;
        for (s = 0; s < *ss; s++) {
            tempsum = 0;
            for (n = 0; n < *nc; n++) {
                tempsum += w[i3(n,s,k,*nc,*ss)];
	    }
            if (tempsum>0) y[k]++;
	}
    }

    *resultcode = 1;  /* generic failure code */
    *value = 0;

    /* long version for now - could easily use frequencies,
       but might later want site-specific and time-specific  */


    if (*type == 0) {
        lambda = *D * M_PI * *sigma * *sigma / 10000;
        for (k=0; k < *kk; k++) {
            tempsum = 0;
            for (j = 0; j <= y[k]; j++) {
                tempsum += choose(y[k],j) * pow(-1, y[k]-j) * 
                    exp(-lambda * (1 - pow(1- *g0, *ss-j) ));
            }
            *value += log(tempsum * choose(*ss, y[k]));
        }
    }
    else {
        par[0] = *g0;
        par[1] = *sigma;
        par[2] = *z;
        for (k = 0; k < *kk; k++) {
            tempsum = 0;
            for (j = 0; j <= y[k]; j++) {
                integral = gintegral2(*fn, *ss - j, par);
                mu = *D * integral / 10000;
                tempsum += choose(y[k],j) * pow(-1, y[k]-j) * exp(-mu);
            }
            *value += log(tempsum * choose(*ss, y[k]));
        }
    }
    *resultcode = 0;
}

