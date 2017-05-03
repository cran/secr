#include "secr.h"

/* unmarked and presence likelihoods */
/* 2011-10-01 */
/* 2011-10-10 pairwise added */

/* temporary unmarked likelihood */
/* only constant model *cc == 1 */
/* likelihood assumes indepndence between sites AND between successive occasions,
   but latter limitation should be fixable */

/*
   can compile with gcc 4.6.3 :
   gcc -Ic:/R/R-3.3.0/include -c unmarked.c -Wall -pedantic -std=gnu99
*/

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
    integral = hintegral(*fn, par);     /* one occasion */
    for (k = 0; k < *kk; k++) {
        for (s = 0; s < *ss; s++) {
            lambda = *D * integral / 10000;
            *value += dpois(nsk[k * *ss + s], lambda, 1); /* log scale */
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
    fnptr fnp = ghnr;
    fnp = getgfnr(fn);
    for (i=0; i<n; i++) {
        x[i] = x[i] * (1 - pow(1 - fnp(tmp,x[i]), tmp[4])); 
    }
}
double hintegral2 (int fn, int sj, double par[]) {
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

/*=======================================================*/
double onepairprb2 (double mua, double mub, double muc,
    double p, int s, int y1, int y2) {
/*=======================================================*/

    /* Calculate probabilities for individual detectors */
    /* Based on Martin Ridout's R code */

    double prob = 0;
    int cmax = 0;
    int c;
    int i;
    double terma, termb;
    double cona, conb;
    double powa, powb;
    double *cprb;

    /* cap cmax to protect against excessive memory demand with rare very large muc */
    cmax = qpois(0.9999, muc, 1, 0);
    if (cmax>1e3)
        cmax = 1e3;
    if (cmax<=0)
        cmax = 0;
    /* Rprintf("cmax %12d  muc %15.8f \n", cmax, muc);  */
    cprb = (double *) R_alloc(cmax+1, sizeof (double));
    for (c=0; c<=cmax; c++)
        cprb[c] = dpois(c, muc, 0);

    for (c=0; c<=cmax; c++) {
        terma = 0;
        for (i=0; i<=y1; i++) {
            cona = choose(y1, i) * pow(-1,y1-i) *
                exp(-mua*(1-pow(1-p,s-i)));
            powa = pow(1-p, s-i);
            terma += cona * pow(powa,c);
        }
        termb = 0;
        for (i=0; i<=y2; i++) {
            conb = choose(y2, i) * pow(-1,y2-i) *
                exp(-mub*(1-pow(1-p,s-i)));
            powb = pow(1-p, s-i);
            termb += conb * pow(powb,c);
        }
        prob += cprb[c] * terma * termb;
    }
    return(choose(s,y1) *  choose(s,y2) * prob);
}
/*==============================================================================*/

void presenceloglik (
    int    *w,           /* capture histories (1:nc, 1:ss, 1:kk) */
    int    *nc,          /* number of rows in w */
    int    *ss,          /* number of occasions */
    int    *kk,          /* number of detectors */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *D,           /* Parameter value - density */
    double *g0,          /* Parameter value - p */
    double *sigma,       /* Parameter value - radius */
    double *z,           /* Parameter value - shape (hazard rate, cumulative gamma etc) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential etc.*/
    int    *type,        /* code 0 = simple, 1 = integrated, 2 = pairwise */
    double *value,       /* return value */
    int    *resultcode   /* 0 if OK */
)

{
    int    i,j,k,n,s;
    int    *y;           /* vector of detector-specific counts, length *kk 0 <= yk <= ss) */
    double d;
    double tempsum;
    double lambda;
    double overlap;
    double sig2;
    double par[3];
    double mu;
    double integral;
    double mua, mub, muc;
 
    /*===============================================================*/

    /*-------------------------------*/
    /* summarise input w as y vector */
    /*-------------------------------*/
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

    /*--------------------------*/
    /* set generic failure code */
    /*--------------------------*/
    *resultcode = 1;
    *value = 0;

    /* long version for now - could easily use frequencies,
       but might later want site-specific and time-specific  */

    if (*type == 0) {

        /*----------------------------------*/
        /* Independent points, fixed radius */
        /*----------------------------------*/

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
    else if (*type == 1) {

        /*------------------------------------------------------*/
        /* Independent points, integrated detection probability */
        /*------------------------------------------------------*/

        par[0] = *g0;
        par[1] = *sigma;
        par[2] = *z;
        for (k = 0; k < *kk; k++) {
            tempsum = 0;
            for (j = 0; j <= y[k]; j++) {
                integral = hintegral2(*fn, *ss - j, par);
                mu = *D * integral / 10000;
                tempsum += choose(y[k],j) * pow(-1, y[k]-j) * exp(-mu);
            }
            *value += log(tempsum * choose(*ss, y[k]));
        }
    }
    else if (*type == 2) {

        /*----------------------------------*/
        /* Pairwise, fixed radius           */
        /*----------------------------------*/

        for (i = 0; i < (*kk - 1); i++) {
            for (j = i+1; j < *kk; j++) {

                /*------------------------------------*/
                /* Distance between detectors i and j */
                /*------------------------------------*/

                d = sqrt( (traps[i] - traps[j]) * (traps[i] - traps[j]) +
                   (traps[i + *kk] - traps[j + *kk]) * (traps[i + *kk] - traps[j + *kk]) );

                /*--------------------------*/
                /* Parameters muA, muB, muC */
                /*--------------------------*/
	        overlap = 0;
                sig2 = *sigma * *sigma;    
                if (d < (2 * *sigma))
                    overlap = 2 * sig2 * acos(d/(2* *sigma)) - d/2 * sqrt(4 * sig2 - d*d);
                mua = (M_PI * sig2 - overlap) * *D / 10000;
                mub = (M_PI * sig2 - overlap) * *D / 10000;
                muc = overlap * *D /10000;

                /*-----------------------------*/
                /* Contribution from this pair */
                /*-----------------------------*/

                *value += log(onepairprb2(mua, mub, muc, *g0, *ss, y[i], y[j]));
            }
            R_CheckUserInterrupt();
        }
    }
    else if (*type == 3) {

        /*--------------------------------*/
        /* Independent points, hazard     */
        /*--------------------------------*/

        par[0] = *g0;
        par[1] = *sigma;
        par[2] = *z;
        for (k = 0; k < *kk; k++) {
            tempsum = 0;
            for (j = 0; j <= y[k]; j++) {
                integral = hintegral2(*fn, *ss - j, par);
                mu = *D * integral / 10000;
                tempsum += choose(y[k],j) * pow(-1, y[k]-j) * exp(-mu);
            }
            *value += log(tempsum * choose(*ss, y[k]));
        }
    }
    else
        error ("unrecognised type");

    *resultcode = 0;
}
/*==============================================================================*/

/* copied from secrdesign */

void Lambda (
    double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 0 = halfnormal */
    double *L,         /* return value vector of length mm */
    int    *resultcode /* 0 for successful completion */
)
{
    int k,m;
    *resultcode = 1;                   /* generic failure */
    if (*fn != 14) error("only hazard halfnormal");
    for (m=0; m<*mm; m++) {
	L[m] = 0;
	for (k = 0; k < *kk; k++) {
	    L[m] += exp(-d2(k, m, traps, mask, *kk, *mm) /2 / par[1] / par[1]);
	}
	L[m] *= par[0];
    }
    *resultcode = 0;                   /* successful completion */
}

void presenceloglik2017 (
    int    *w,           /* capture histories (1:nc, 1:ss, 1:kk) */
    int    *nc,          /* number of rows in w */
    int    *ss,          /* number of occasions */
    int    *kk,          /* number of detectors */
    int    *mm,          /* number of mask points */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *mask,        /* x,y locations of mask points (first x, then y) */
    double *cellarea,
    double *D,           /* Parameter value - density */
    double *lambda0,     /* Parameter value - lambda0 */
    double *sigma,       /* Parameter value - radius */
    double *z,           /* Parameter value - shape (hazard rate, cumulative gamma etc) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential etc.*/
    double *value,       /* return value */
    int    *resultcode   /* 0 if OK */
)

{
    int    k,n,s;
    int    *y;           /* vector of detector-specific counts, length *kk 0 <= yk <= ss) */
    double *kappa;       /* vector of detector-specific kappa, length *kk */
    double tempsum;
    double p;
    double par[3];
 
    /*===============================================================*/

    /*-------------------------------*/
    /* summarise input w as y vector */
    /*-------------------------------*/
    y = (int *) R_alloc(*kk, sizeof(int));
    kappa = (double *) R_alloc(*kk, sizeof(double));
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

    /*--------------------------*/
    /* set generic failure code */
    /*--------------------------*/
    *resultcode = 1;
    *value = 0;
    par[0] = *lambda0;
    par[1] = *sigma;
    par[2] = *z;

    /*-----------------------------------*/
    /* Independent points, summed hazard */
    /*-----------------------------------*/

    Lambda(par, mm, kk, mask, traps, fn, kappa, resultcode);
    for (k=0; k < *kk; k++) {
	p = 1 - exp(-kappa[k] * *ss * *D * *cellarea);
	if (y[k]>0)
	    *value += log (p);
	else
	    *value += log (1-p);
    }

    *resultcode = 0;
}
/*==============================================================================*/

