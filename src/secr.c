/*
   External procedures for secr package

   can compile with gcc 4.2.1-sjlj :
   gcc -Ic:/R/R-2.13.0/include -c secr.c -Wall -pedantic -std=gnu99

   [confirmed 2011-06-13 0059]

*/
/* 2011-03-19 */
/* pmix does not work with nc = 0 */

/* 2011-04-05 allow 'partial likelihood' option in secrloglik */
/* *like == 2 */
/* 2011-05-05 tweak integralprw1 to admit factor of D */
/* integralprw1 acquires *useD argument 0/1 */
/* D may be provided as a third column of *mask */
/* 2011-05-12 clean l 726, 4136 */
/* 2011-06-07 alongtransect */
/* 2011-06-21 SegCircle2; debugged 06-22 */
/* 2011-06-21 mod integral1D and gintegral1 to handle fn == 4 differently */
/* 2011-06-22 changed all bitwise logical operators to &&, || */
/* 2011-06-22 changed all floating point abs to fabs */
/* 2011-09-30 major reorgansiation spinning off other '.c' files */
/* 2011-09-30 re-write of pfn related to above*/
/* 2011-11-15 parameter index array renamed 'gsb' to PIA' for clarity */
/* 2011-11-15 substantial revision of exclusive detector code */
/* 2012-01-22 removed turnover (phi) code Only in 2.3.1 */
/* 2012-01-30 prwisignal modified to allow missing signal strengths (<=0) */
/* 2012-01-31 cuerate replaced by miscparm, which may be a vector of parameters */
/* 2012-02-01 also collaps mufn,mufnsph to one and simplified dnorm call */

#include "secr.h"
#include <time.h>

FILE *out;      /* for debugging */
/*==============================================================================*/

void R_CheckUserInterrupt(void);

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

double gbinom(int count, int size, double p, int uselog)
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
    else
        return (dbinom (count, size, p, uselog));
}

double gnbinom (int count, int size, double mu, int uselog)
{

    /* prob = size / (size + mu) */
    size = fabs(size);  /* in case negative 'binomN' passed */

    if (count == 0) {  /* faster - added 2010-10-11*/
        if (uselog) return( log(size/(size+mu)) * log(size) );
        else return (pow(size/(size+mu), size));
    }
    else
        return (dnbinom (count, size, size/(size+mu), uselog));
}

double gbinomFP (int count, double size, double p, int uselog)
{
    return ( lgamma(size+1) - lgamma(size-count+1) - lgamma(count+1) +
             count * log(p) + (size - count) * log (1-p) );
}

/*==============================================================================*/

/*===============================================================*/
void gxy (int *n, int *fn, double *par, double *w, double *xy) {
    int maxj = 1000000;
    double r;
    double theta;
    fnptr fnp = hn;
    int i = 0;
    int j;
    fnp = gethfn(*fn);
    for (i=0; i< *n; i++) {
        theta = unif_rand() * 2 * M_PI;
        for (j=0; j<maxj; j++) {
            r = *w * sqrt(unif_rand());
            if (unif_rand() < fnp(par, r))
                break;
        }
        xy[i]      = r * cos(theta);
        xy[*n + i] = r * sin(theta);
    }
}

/* find upper and lower points on perimeter of poly at x-coordinate x */
void yab(double x[], int *i, int *np, double poly[], double *a, double *b) {
    int k;
    int nv = 0;
    double ab[3];
    /* note 'sign' is RMath function */
    for (k=0; k< (*np-1); k++) {
        if (sign(poly[k]- x[*i]) != sign(poly[k+1]- x[*i])) {
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
    fnptr fnp = hn;
    p = (double*) ex;
    fn = round(p[3]);
    mx = p[4];
    my = p[5];
    xy[0] = p[6];

    /* set detection function */
    fnp = gethfn(fn);
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
    int reportier = 0;

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
    if ((ier != 0) & (reportier))
        Rprintf("ier error code in integral2D %5d\n", ier);
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

    /* limits from bounding box of this polygon */
    ns = *n2 - *n1 + 1;
    for (k=0; k<ns; k++) {
        ax = fmin2(ax, traps[k+ns]);
        bx = fmax2(bx, traps[k+ns]);
    }
    ex = (double *) R_alloc(10 + 2 * *kk, sizeof(double));
    ex[0] = gsbval[*c];   /* 1.0? */
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

/* probability of count */

double countp (int count, int binomN, double lambda) {
    /* Poisson */
    if (binomN == 0)
        return ( gpois (count, lambda, 0));

    /* Bernoulli */
    else if (binomN == 1)
        if (count == 0)
            return ( 1 - lambda );
        else
            return ( lambda );

    /* negative binomial */
    else if (binomN < 0)
        return ( gnbinom (count, binomN, lambda, 0) );

    /* binomial */
    else
        return ( gbinom (count, binomN, lambda / binomN, 0) );
}
/*=============================================================*/

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
    fnptr fnp = hn;
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
    fnp = gethfn(fn);
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

    /* 2011-06-21 uniform treated separately */
    if (fn == 4) {
        for (k=n1+1; k<=n2; k++) {  /* upper bound is length of this transect */
            bx += SegCircle2(
               traps[k-1], traps[k-1+kk],
               traps[k], traps[k+kk],
               mask[m], mask[m+mm],
	       gsbval[cc + c]);
        }
        return (bx);
    }

    for (k=n1+1; k<=n2; k++) {  /* upper bound is length of this transect */
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

double pndot (int m, int n, int s1, int s2, int x, int ncol, int PIA0[],
	      double gk0[], int ss, int kk, int cc0, int nmix,
              double gsb0val[], int param)
/*
    probability animal at point m on mask is caught
    n may indicate group (full likelihood; ncol= number of groups) or
    individual (conditional likelihood; ncol= number of individuals)
    aligned with secrloglik 2009 06 25

    2009 10 24 adjusted to allow summation over qq < ss
    2009 11 12 'kk' should be number of parts for polygon detectors
    2011 01 04 'param' for GR parameterisation of multi
*/
{
    int k,s,c;
    int gi;
    int wxi;
    double pp;
    double p0;
    double Ei;
    pp = 1;
    if (param == 0) {
        for (s=s1-1; s<s2; s++) {
            for (k=0; k< kk; k++) {
                wxi = i4(n,s,k,x,ncol,ss,kk);
                c = PIA0[wxi] - 1;
                if (c >= 0) {    /* drops unset traps */
                    gi = i3(c,k,m,cc0,kk);
                    pp *= 1 - gk0[gi];
                }
            }
        }
        return (1 - pp);
    }
    else {
        for (s=s1-1; s<s2; s++) {
            Ei = 0;
            wxi = i4(n,s,0,x,ncol,ss,kk);
            c = PIA0[wxi] - 1;
            p0 = gsb0val[c];     /* value for first detector */
            for (k=0; k< kk; k++) {
                wxi = i4(n,s,k,x,ncol,ss,kk);
                c = PIA0[wxi] - 1;
                if (c >= 0) {    /* drops unset traps */
                    if (p0 != gsb0val[c])
			error("trap-specific p0 not allowed with G&R param");
                    gi = i3(c,k,m,cc0,kk);
                    Ei += gk0[gi] / p0;
                }
            }
            pp *= (1 - p0 * exp(-1/Ei));
        }
        return (1 - pp);
    }
}
/*===============================================================*/

void integralprw1 (
    int    *detect,    /* detector 0 multi, 1 proximity etc. */
    int    *param,     /* parameterisation 0 Borchers & Efford 1 Gardner & Royle */
    double *gsb0val,   /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *nc,        /* number of individuals */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    int    *nmix,      /* number of mixtures */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *cc0,       /* number of g0/sigma/b combinations [naive animal] */
    int    *PIA0,      /* lookup which g0/sigma/b combination to use for given n, S, K
                          [naive animal] */
    int    *ncol,      /* number of columns in PIA0; added 2009 06 25 */
    double *area,      /* area associated with each mask point (ha) */
    double *miscparm,  /* miscellaneous parameters (cuerate etc) */
    int    *fn,        /* codes
                          0 = halfnormal,
                          1 = hazard rate,
                          2 = exponential,
                          9 = binary signal strength,
                          10 = signal strength,
                          11 = signal strength spher, */
    int    *binomN,    /* number of trials for 'count' detector modelled with binomial */
    int    *useD,      /* logical : does third column of mask contain D weight? 2011-05-05 */
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
    double hk;
    double par0[4];
    double D = 1.0;
    int cumk[1001];
    int nk = 0;
    double *pmix;      /* proportion in each mixture */
    int gpar = 2;
    double stdint;

    *resultcode = 1;                   /* generic failure */

    /* determine number of polygons or transects */
    /* for polygon or transects detectors, kk is vector ending in zero */
    if ((*detect==3) || (*detect==4) || (*detect==6) || (*detect==7)) {
        cumk[0] = 0;
        for (i=0; i < maxnpoly; i++) {
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
    }
    else nk = *kk;

    if ((*detect==0) || (*detect==1) || (*detect==3) || (*detect==4))
        *binomN = 1;

    /* Allocate space for array of naive detection probability */
    gk0 = (double *) R_alloc(*cc0 * nk * *mm, sizeof (double));
    pmix = (double *)  R_alloc(*nc * *nmix, sizeof (double));
    for (i=0; i< *nc * *nmix; i++) pmix[i] = 1; /* default */

    gfn = getgfn(*fn);

    if (*nmix>1) gpar++;

    if (*detect == 0) {
        for (c=0; c<*cc0; c++) {
            for (k=0; k<*kk; k++) {
                for (m=0; m<*mm; m++) {
                        gi = i3(c, k, m, *cc0, nk);
                        gk0[gi] = gfn(k, m, c, gsb0val, *cc0, traps, mask, nk, *mm, miscparm);
                }
            }
        }
    }
    else if ((*detect == 1) || (*detect == 2) || (*detect == 5) || (*detect==8)|| (*detect==9)) {
        for (c=0; c<*cc0; c++) {
            for (k=0; k<*kk; k++) {
                for (m=0; m<*mm; m++) {
                    lambda = gfn(k, m, c, gsb0val, *cc0, traps, mask, nk, *mm, miscparm);
                    gi = i3(c, k, m, *cc0, nk);
                    gk0[gi] = 1 - countp (0, *binomN, lambda);
                }
            }
        }
    }
    else if ((*detect == 3) || (*detect == 6)) {
        for (c=0; c<*cc0; c++) {
            par0[0] =  gsb0val[c];
            par0[1] = gsb0val[*cc0 + c];
            par0[2] = gsb0val[2* *cc0 + c];
            stdint = gintegral(*fn, par0);   /* DOES THIS ASSUME NO TRAP(POLYGON) EFFECT ? */
            for (k=0; k<nk; k++) {                 /* over polygons */
                for (m=0; m<*mm; m++) {
                    hk = par0[0] * integral2D (*fn, m, c, gsb0val, *cc0,
                        traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
                    gi = i3(c, k, m, *cc0, nk);
                    gk0[gi] = 1 - countp (0, *binomN, hk);
                }
            }
        }
    }
    else if ((*detect == 4) || (*detect == 7)) {
        for (c=0; c<*cc0; c++) {
            par0[0] = gsb0val[c];
            par0[1] = gsb0val[*cc0 + c];
            par0[2] = gsb0val[2* *cc0 + c];
            stdint = gintegral1(*fn, par0);   /* DOES THIS ASSUME NO TRAP(POLYGON) EFFECT ? */
            for (k=0; k<nk; k++) {                 /* over transects */
                for (m=0; m<*mm; m++) {
                    hk = par0[0] * integral1D (*fn, m, c, gsb0val, *cc0, traps,
                        mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
                    gi = i3(c, k, m, *cc0, nk);
                    gk0[gi] = 1 - countp (0, *binomN, hk);
                }
            }
        }
    }
    else error ("unrecognised detector type in external C fn integralprw1");

    if (*nmix>1) {
        for (n=0; n<*nc; n++) {
            for (x=0; x<*nmix; x++) {
                wxi = i4(n,0,0,x,*ncol,*ss,*kk);
                c = PIA0[wxi] - 1;
                pmix[*nmix * n + x] = gsb0val[*cc0 * (gpar-1) + c];
           }
       }
    }

    /* in case of zero detections 2011-03-19 */
    if (*nc == 0)
        *nc = 1;

    for (n=0; n<*nc; n++) {            /* CH numbered 0 <= n < *nc */
        if ((*ncol > 1) || (n == 0)) {  /* no need to repeat if constant */
            if ((n+1) > *ncol) {       /* groups not implemented */
                *resultcode = 3;
                return;
            }
            asum = 0;
            for (x=0; x<*nmix; x++) {
                for (m=0; m<*mm; m++) {
                    if (*useD>0) D = mask[2 * *mm + m];
                    asum += pmix[*nmix * n + x] * D * pndot (m, n, 1, *ss, x,
                        *ncol, PIA0, gk0, *ss, nk, *cc0, *nmix, gsb0val, *param);
                }
            }
        }
        a[n] = *area * asum;
    }
    *resultcode = 0;                   /* successful completion */
}

/*==============================================================================*/

double prwimulti

   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc,
    int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], 
    double mask[], double minp)

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
    int debug = 0;
    for (s=s1-1; s<s2; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        if (htemp < fuzz) { result = 0; break; }
        k = w[nc * s + n];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        if (k > 0) {
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = -log (1 - gk[gi]);
            pks *= (1-expmin(-htemp)) / htemp;
	    if((m==1311) && debug)
	    Rprintf("n %5d  k %5d  c %5d, gi %5d, pks %12.8f \n", n, k, c, gi, pks); 
        }
        else  {
            pks = expmin(-htemp);      /* Not captured */
        }
         result *= pks;
        if (dead) break;
    }
    return (result);
}
/*=============================================================*/

double prwimultiGR
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk, int ss,
    int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < *nc)
    given that animal's range centre is at m
    MULTI-CATCH DETECTOR, GARDNER & ROYLE PARAMETERISATION
*/
{
    int s;                             /* index of occasion  0 <= s < *ss */
    int k;                             /* index of trap      0 <= k < *kk */
    int c;
    int gi;
    int dead = 0;
    double pks;
    double result = 1.0;
    double htemp;
    double Ei;
    double p0;

    for (s=s1-1; s<s2; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        k = w[nc * s + n];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        if (k > 0) {
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            p0 = gsbval[c];
            Ei = htemp / p0;
            if (Ei < fuzz) { result = 0; break; }
            gi = i3(c, k-1, m, cc, kk);
            pks = p0 * exp(-1/Ei);      /* captured */
            pks *= gk[gi]/p0 / Ei;      /* captured in trap k */
        }
        else  {
            c = PIA[i4(n, s, 0, x, nc, ss, kk)] - 1;   /* use combination for k=0 (good as any) */
            p0 = gsbval[c];
            Ei = htemp / p0;
            if (Ei < fuzz) { result = 0; break; }
            pks = 1 -  p0 * exp(-1/Ei);      /* Not captured */
        }
        result *= pks;
        if (dead) break;
    }
    return (result);
}
/*=============================================================*/

double prwicount
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
     double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk,
     int ss, int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], 
     double minp)

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
            c = PIA[wxi] - 1;
            if (c >= 0) {                                       /* skip if this trap not set */
                gi  = i3(c, k, m, cc, kk);
                result *= countp (count, binomN, gk[gi]);
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
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk,
    int ss, int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[],
    double minp)

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    SIGNAL STRENGTH DETECTOR
    Modified 2012-01-30 to allow missing signal strengths
    Does not need cut; simplified dnorm call 2012-02-01
*/

{
    int s;   /* index of occasion  0 <= s < ss */
    int k;   /* index of trap      0 <= k < kk */
    int j;   /* index of detection */
    int c, wi, wxi, gi;
    double result = 1.0;
    double mu, sdS;
    int start = 0;
    int count = 0;
    int spherical;
    double sig = 0;

    spherical = round (detspec[3]);
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<kk; k++) {
            wxi = i4(n, s, k, x, nc, ss, kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                                           /* skip if this trap not set */
                wi = i3(n, s, k, nc, ss);
                count = abs(w[wi]);                                 /* ignore 'deaths */
                if (count==0) {
                    /* not detected at this mic */
                    gi = i3(c,k,m,cc,kk);
                    result *= countp (0, binomN, gk[gi]);
                }
                else {
                    /* detected at this mic */
                    start = round(detspec[wi+4]);
                    for (j=0; j<count; j++) {
			sig = signal[start+j];
			if (sig >= 0) {
                            /* valid measurement of signal */
			    mu = mufn (k, m, gsbval[c], gsbval[cc + c],
				       traps, mask, kk, mm, spherical);
			    sdS = gsbval[cc * 2 + c];
			    result *= dnorm((sig - mu), 0, sdS, 0);
			}
			else  {
                            /* signal value missing; detection only */
			    gi = i3(c,k,m,cc,kk);
			    result *= countp (1, binomN, gk[gi]);
			}			   
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

double prwisignalnoise
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk,
    int ss, int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[],
    double minp)

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    SIGNAL + NOISE DETECTOR 2012-02-07

detspec[0] nk
detspec[1] nd
detspec[2] miscparm[0]  cut
detspec[3] miscparm[1]  muN
detspec[4] miscparm[2]  sdN
detspec[5] spherical
detspec[6] start[0] (i=0, s=0, k=0)
..
detspec[6+nc*ss*nk-1] start[nc*ss*nk-1]

start[w_isk] is index to signal for w_isk 
start[w_isk] + nd is index to noise for w_isk 

*/

{
    int s;   /* index of occasion  0 <= s < ss */
    int k;   /* index of trap      0 <= k < kk */
    int j;   /* index of detection */
    int c, wi, wxi, gi;
    double result = 1.0;
    int start = 0;
    int count = 0;
    int spherical;
    double sig = 0;
    double nois = 0;
    double muS;
    double sdS;
    double muN;
    double sdN;
    /* double cut; */
    int nd;
    /* double xi; */
    double f = 1;

    /* xi = 10 /  M_LN10; */
    nd = round (detspec[1]);
    spherical = round (detspec[5]);
    /* cut = detspec[2]; */
    muN = detspec[3];
    sdN = detspec[4];
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<kk; k++) {
            wxi = i4(n, s, k, x, nc, ss, kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                                           /* skip if this trap not set */
                wi = i3(n, s, k, nc, ss);
                count = abs(w[wi]);                                 /* ignore 'deaths */
                if (count==0) {
                    /* not detected at this mic */
                    gi = i3(c,k,m,cc,kk);
                    result *= countp (0, binomN, gk[gi]);
                }
                else {
                    /* detected at this mic */
                    start = round(detspec[wi+6]);
                    for (j=0; j<count; j++) {
			sig = signal[start+j];
			nois = signal[start+nd+j];

			if (sig >= 0) {
                            /* valid measurement of signal */
			    muS = mufn (k, m, gsbval[c], gsbval[cc + c],
				       traps, mask, kk, mm, spherical);
                            /* acknowledge noise component of observed signal: SLOW
                               muS = xi * log(pow(10, muS/10) + pow(10, muN/10));  */
			    sdS = gsbval[cc * 2 + c];
			    f =  dnorm(sig - nois, muS-muN, sqrt(sdS*sdS+sdN*sdN), 0) *
				dnorm (nois, muN, sdN,0);   /* does not allow for effect of    */
                                                            /* truncation on N distribution !! */
 			    result *= f; 
			}
			else  {
                            /* signal value missing; detection only */
			    gi = i3(c,k,m,cc,kk);
			    result *= countp (1, binomN, gk[gi]);
			}			   
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

void fnst (double *x, int n, void *ex){
    int i;
    double * st;
    st = (double*) ex;   /* [s, t] */
    for(i=0; i<n; i++) 
	x[i] = exp(st[0] * x[i] - x[i]*x[i]/2) * pow(pnorm(x[i],0,1,1,0), st[1]-1);
}

double d2lnorm (double x, double mu1, double s1, double mu2, double s2) {
/* approximate density of sum of two lognormal random variables */
/* Szyszkowiccz and Yanikomeroglu */
    double s;
    double m;
    double t;
    double y;
    double Lambda = 0;
    double ex[2];
    double bound = 0;
    int    inf = 2;
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

    s = fmax2(s1,s2);
    t = s*s * (1/(s1*s1) + 1/(s2*s2));
    if (s1 == s2) {
        m = log(exp(mu1) + exp(mu2)) - log(2) - pnorm(s/sqrt(2),0,1,1,1);
    }
    else {
        ex[0] = s;
        ex[1] = t;
        Rdqagi(fnst, ex, &bound, &inf, &epsabs, &epsrel, &Lambda, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
	m = log (exp(mu1 + s1*s1/2) + exp(mu2 + s2*s2/2)) - log(Lambda) - log(t) + log(M_2PI)/2;
    }
    y = (log(x) - m) / s;
    return (t / (sqrt(M_2PI) * x * s) * exp( - y * y / 2 ) * pow( pnorm (y,0,1,0,0), t-1));
}

double prwisignal2
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk,
    int ss, int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[],
    double minp)

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    SIGNAL STRENGTH DETECTOR
    Modified 2012-01-30 to allow missing signal strengths
    Does not need cut: simplified dnorm call 2012-02-01

    Model sum of two lognormals
*/

{
    int s;   /* index of occasion  0 <= s < ss */
    int k;   /* index of trap      0 <= k < kk */
    int j;   /* index of detection */
    int c, wi, wxi, gi;
    double result = 1.0;
    double mu, sdS;
    int start = 0;
    int count = 0;
    int spherical;
    double sig = 0;
    double xi;
    double f = 1;

    xi = 10 / log(10);

    spherical = round (detspec[3]);
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<kk; k++) {
            wxi = i4(n, s, k, x, nc, ss, kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                                           /* skip if this trap not set */
                wi = i3(n, s, k, nc, ss);
                count = abs(w[wi]);                                 /* ignore 'deaths */
                if (count==0) {
                    /* not detected at this mic */
                    gi = i3(c,k,m,cc,kk);
                    result *= countp (0, binomN, gk[gi]);
                }
                else {
                    /* detected at this mic */
                    start = round(detspec[wi+4]);
                    for (j=0; j<count; j++) {
			sig = signal[start+j];
			if (sig >= 0) {
                            /* valid measurement of signal */
			    mu = mufn (k, m, gsbval[c], gsbval[cc + c],
				       traps, mask, kk, mm, spherical);
			    sdS = gsbval[cc * 2 + c];
                            /* fixed mu2, s2 : ovenbird noise */
			    f = d2lnorm(sig/xi, mu/xi, sdS/xi, detspec[1]/xi, detspec[2]/xi);
			    if (m == 2200)
/*	    Rprintf("sig %12.6f mu %12.6f sdS %12.6f f %12.6f \n", sig,mu,sdS,f); */
 			    result *= f; 
			}
			else  {
                            /* signal value missing; detection only */
			    gi = i3(c,k,m,cc,kk);
			    result *= countp (1, binomN, gk[gi]);
			}			   
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
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk,
    int ss, int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[],
    double minp)
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
            c = PIA[wxi] - 1;
            if (c >= 0) {                       /* skip if this trap not set */
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
                result *= countp (0, binomN, lambda);
                if (count>0) {
                    start = round(detspec[wi]);  /* plus offset ? */
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
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk, int ss,
    int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)

/*
    Likelihood component due to capture history n (0 <= n < nc)
    given that animal's range centre is at m
    POLYGON DETECTOR
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of part 0 <= k < nk  */
    int j;   /* index of xy record */
    int c, wxi, wi, gi;
    long count;
    int dead = 0;
    double result = 1.0;
    int nk;
    int nd;
    int start = 0;
    double H;

    nk = round(detspec[0]);
    nd = round(detspec[1]);
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<nk; k++) {
            wi = i3(n,s,k,nc,ss);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            wxi = i4(n,s,k,x,nc,ss,kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                          /* skip if this polygon not used */
                gi  = i3(c,k,m,cc,nk);
                result *= countp (count, binomN, gk[gi]);
                if ((result > minp) && (count>0)) {               /* avoid underflow 2009 11 16 */
                    start = round(detspec[wi+cc+2]);
                    for (j=start; j < start+count; j++) {
                        /* convoluted way of retrieving integral2D(gfn(x) over k)) */
			H = gk[gi] / gsbval[c] * detspec[2+c] / 10000;  /* why 10000? */
                        result *= gfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / H;
                    }
		}
            }
        }
        if (result <= minp) {
            /* Pr(detection history) fell below minprob in prwipolygon */
            /* Simply aborting at this point does not work 2011-01-30 */
            /* result = -1; break; */
            result = minp; break;
        }
        if (dead==1) break;
    }
    return (result);
}
/*=============================================================*/

double prwipolygonX
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk, int ss,
    int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < *nc)
    given that animal's range centre is at m
    EXCLUSIVE POLYGON DETECTOR
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
    int j, wi;
    int nd;
    double H;

    nd = round(detspec[1]);

    for (s=s1-1; s<s2; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        if (htemp < fuzz) { result = 0; break; }
        wi = nc * s + n;
        k = w[wi];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        /* PR | detected */
        if (k > 0) {
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = -log (1 - gk[gi]);
            pks *= (1-expmin(-htemp)) / htemp;
        }
        /* PR | not detected */
        else  {
            pks = expmin(-htemp);      /* Not captured */
        }
        result *= pks;
        /* Pr(location) */
        if (k > 0) {
            j = round(detspec[wi+cc+2]);
            H = gk[gi] / gsbval[c] * detspec[2+c] / 10000;
            result *= gfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / H;
        }

        if (dead) break;
    }
    return (result);
}
/*=============================================================*/

double prwitransectX
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk, int ss,
    int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < *nc)
    given that animal's range centre is at m
    EXCLUSIVE TRANSECT DETECTOR
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
    int j, wi;
    int nd;
    double H;

    nd = round(detspec[1]);

    for (s=s1-1; s<s2; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        if (htemp < fuzz) { result = 0; break; }
        wi = nc * s + n;
        k = w[wi];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        /* PR | detected */
        if (k > 0) {
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = -log (1 - gk[gi]);
            pks *= (1-expmin(-htemp)) / htemp;
        }
        /* PR | not detected */
        else  {
            pks = expmin(-htemp);      /* Not captured */
        }
        result *= pks;
        /* Pr(location) */
        if (k > 0) {
            j = round(detspec[wi+cc+2]);
            H = gk[gi] / gsbval[c] * detspec[2+c] / 100;
            result *= gfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / H;
        }

        if (dead) break;
    }
    return (result);
}
/*=============================================================*/

double prwitransect
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int PIA[],
    double gk[], int binomN, double detspec[], double h[], int hindex[], int cc, int nc, int kk, int ss,
    int mm, int nmix, gfnptr gfn, double gsbval[], double traps[], double mask[], double minp)

/*
    Likelihood component due to capture history n (0 <= n < nc)
    given that animal's range centre is at m
    TRANSECT DETECTOR
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of part 0 <= k < nk  */
    int j;   /* index of xy record */
    int c, wxi, wi, gi;
    long count;
    int dead = 0;
    double result = 1.0;
    int nk;
    int nd;
    int start = 0;
    double H;

    nk = round(detspec[0]);
    nd = round(detspec[1]);
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<nk; k++) {
            wi = i3(n,s,k,nc,ss);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            wxi = i4(n,s,k,x,nc,ss,kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                          /* skip if this transect not used */
                gi  = i3(c,k,m,cc,nk);
                H = gk[gi];
                result *= countp (count, binomN, H);
                if ((result > minp) && (count>0)) {               /* avoid underflow 2009 11 16 */
                    start = round(detspec[wi+cc+2]);
                    for (j=start; j < start+count; j++) {
                        /* convoluted way of retrieving integral1D(gfn(x) over k)) */
			H = gk[gi] / gsbval[c] * detspec[2+c] / 100;   /* why 100? */
                        result *= gfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / H;
                    }
		}
            }
        }
        if (result <= minp) {
            /* Pr(detection history) fell below minprob in prwipolygon */
            /* Simply aborting at this point does not work 2011-01-30 */
            /* result = -1; break; */
            result = minp; break;
        }
        if (dead==1) break;
    }
    return (result);
}
/*=============================================================*/

void pdotpoint (double *xy, int *nxy, double *traps, int *used, int *kk,
		int *fn, double *par, int *nocc, double *w2, int *binomN,
                double *value)
{
    int i;
    int k;
    double dk2;
    double tempval;
    double g0;
    double sigma;
    double z = 1;
    double cutval [1];
    int s;
    int sumused = 0;
    cutval[0] = 0;
    for (k=0; k<*kk; k++) {
        for (s=0; s<*nocc; s++) {
            sumused += used[s * *kk + k];
        }
    }

    if (*fn>13)
        if (*fn != 20) error("require detectfn<14");
    g0 = par[0];
    sigma = par[1];
    if (!((*fn == 0) || (*fn == 2) || (*fn == 4) || (*fn == 9) || (*fn == 20)))
        z = par[2];
    if ((*fn == 10) || (*fn == 11) || (*fn == 12) || (*fn == 13))
        cutval[0] = par[3];

    if (sumused < (*nocc * *kk)) {
        for (i=0; i<*nxy; i++) {
            tempval = 1;
            for (s=0; s<*nocc; s++) {
                for (k=0; k<*kk; k++) {
                    if (used[s * *kk + k] > 0) {
                        dk2 = (xy[i]-traps[k]) * (xy[i]-traps[k]) +
                              (xy[i + *nxy]-traps[k+ *kk]) * (xy[i + *nxy]-traps[k+ *kk]);
                        tempval *= 1 - pfn(*fn, dk2, g0, sigma, z, cutval, *w2);
                    }
                }
            }
            value[i] = 1 - tempval;
        }
    }
    else {
        if (*fn == 20) {    /* Gardner et al 2009 */
            for (i=0; i<*nxy; i++) {
                tempval = 0;
                for (k=0; k<*kk; k++) {
                    dk2 = (xy[i]-traps[k]) * (xy[i]-traps[k]) +
                          (xy[i + *nxy]-traps[k+ *kk]) * (xy[i + *nxy]-traps[k+ *kk]);
                    tempval += pfn(0, dk2, 1, sigma, z, cutval, *w2);
                }
                if (tempval>0)
                    value[i] = 1 - pow (1 - g0 * exp(-1/tempval), *nocc);
                else
                    value[i] = 0;
            }
        }
        else {
            for (i=0; i<*nxy; i++) {
                tempval = 1;
                for (k=0; k<*kk; k++) {
                    dk2 = (xy[i]-traps[k]) * (xy[i]-traps[k]) +
                          (xy[i + *nxy]-traps[k+ *kk]) * (xy[i + *nxy]-traps[k+ *kk]);
                    tempval *= 1 - pfn(*fn, dk2, g0, sigma, z, cutval, *w2);
                }
                value[i] += 1 - pow(tempval, *nocc);
            }
        }
    }
}

/*=============================================================*/

void pdotpoly (double *xy, int *nxy, double *traps, int *used, int *nk,
    int *kk, int *fn, double *par, int *nocc, int *binomN, double *value)
{
    int k;
    int i;
    double hk = 0;
    double sumhk;
    double tempval;
    int *cumk;
    int s;
    int sumused = 0;
    double stdint;

    for (k=0; k<*nk; k++) {
        for (s=0; s<*nocc; s++) {
            sumused += used[s * *nk + k];
        }
    }
    cumk = (int *) R_alloc(*nk+1, sizeof(int));
    cumk[0] = 0;
    for (i=0; i<*nk; i++)
        cumk[i+1] = cumk[i] + kk[i];

    stdint = gintegral(*fn, par);

    /* was *kk 2011-01-24 */
    if (sumused < (*nocc * *nk)) {
        for (i=0; i<*nxy; i++) {
            tempval = 1.0;
            for (s=0; s<*nocc; s++) {
                sumhk = 0.0;
                for (k=0; k<*nk; k++) {
                    if (used[s * *nk + k] > 0) {
                        hk = par[0] * integral2D (*fn, i, 0, par, 1, traps, xy,
                            cumk[k], cumk[k+1]-1, cumk[*nk], *nxy) / stdint;
                        sumhk += hk;
                    }
                }
                tempval *= countp(0, *binomN, sumhk);
            }
            value[i] = 1 - tempval;
        }
    }
    else {
        for (i=0; i<*nxy; i++) {
            sumhk = 0.0;
            for (k=0; k<*nk; k++) {                        /* over parts */
                hk = par[0] * integral2D (*fn, i, 0, par, 1, traps, xy,
                    cumk[k], cumk[k+1]-1, cumk[*nk], *nxy) / stdint;
                sumhk += hk;
            }
            tempval = countp(0, *binomN, sumhk);
            value[i] = 1 - pow(tempval, *nocc);
        }
    }
}
/*=============================================================*/

void pdottransect (double *xy, int *nxy, double *traps, int *used, int *nk,
    int *kk, int *fn, double *par, int *nocc, int *binomN, double *value)
{
    int k;
    int i;
    double hk;
    double sumhk;
    double tempval;
    int *cumk;
    int s;
    int sumused = 0;
    double stdint;

    for (k=0; k<*nk; k++) {
        for (s=0; s<*nocc; s++) {
            sumused += used[s * *nk + k];
        }
    }
    cumk = (int *) R_alloc(*nk+1, sizeof(int));
    cumk[0] = 0;
    for (i=0; i<*nk; i++)
        cumk[i+1] = cumk[i] + kk[i];

    stdint = gintegral1(*fn, par);

    if (sumused < (*nocc * *nk)) {
        for (i=0; i<*nxy; i++) {
            tempval = 1.0;
            for (s=0; s<*nocc; s++) {
                sumhk = 0.0;
                for (k=0; k<*nk; k++) {                        /* over transects */
                    if (used[s * *nk + k] > 0) {
                        hk = par[0] * integral1D (*fn, i, 0, par, 1, traps, xy,
                            cumk[k], cumk[k+1]-1, cumk[*nk], *nxy) / stdint;
                        sumhk += hk;
                    }
                }
                tempval *= countp(0, *binomN, sumhk);
            }
            value[i] = 1 - tempval;
        }
    }
    else {
        for (i=0; i<*nxy; i++) {
            sumhk = 0.0;
            for (k=0; k<*nk; k++) {                        /* over transects */
                hk = par[0] * integral1D (*fn, i, 0, par, 1, traps, xy,
                    cumk[k], cumk[k+1]-1, cumk[*nk], *nxy) / stdint;
                sumhk += hk;
            }
            tempval = countp(0, *binomN, sumhk);
            value[i] = 1 - pow(tempval, *nocc);
        }
    }
}

/*=============================================================*/
/* 2011-01-11 */
double pndotgrp (double pd, double cuerate, double tol)
/* pd = Pr(cue detected) */
{
    int maxcount = 100;
    int count = 1;
    double dp = 1;
    double psum = 0;
    if (tol <= 0)
        error("requires positive tol");
    count = 1;
    while ((count < maxcount) && (dp > tol)) {
        dp = dpois(count, cuerate, 0);
        psum += dp * (1 - pow(pd, count));
        count++;
    }
    return (1 - psum);
}
/*=============================================================*/

double hazard (double pp) {
    if (pp > (1-fuzz))  /* pp close to 1.0 - approx limit */
        pp = huge;      /* g0 very large (effecti inf hazard) */
    else {
        if (pp <= 0) 
            pp = 0;
        else 
            pp = -log(1-pp);
    }
    return(pp);
}
/*=============================================================*/

void countsessions(int *jj, int J[], int ss, double intervals[]) {
    int s;

    /* index of primary session corresp secondary session s */
    J[0] = 0;
    for (s = 1; s < ss; s++) {
        if (intervals[s-1] > 1e-10)  {
            J[s] = J[s-1] + 1;
        }
        else {
            J[s] = J[s-1];
        }
    }
    /* number of primary sessions */
    *jj = J[ss-1] + 1;
}
/*=============================================================*/

void secrloglik (
    int    *like,        /* likelihood 0 full, 1 conditional */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
    int    *param,       /* parameterisation 0 Borchers & Efford 1 Gardner & Royle */
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
    int    *PIA,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *PIA0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive] */
    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameters (cuerate, cutval etc.) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *minprob,     /* minimum value of P(detection history) */
    double *a,           /* a(theta) */
    double *value,       /* return value integral */
    int    *resultcode   /* 0 if OK */
)

/*
  A note on ordering of input data arrays 2012-02-09

  Using 'i' to subscript animals, 's' to subscript occasions, 'k' to subscript detectors..

 'w' is in isk order - dim(CH) = c(nc,ss,nk) or c(nc,ss) for 2-D exclusive detectors 

 'signal' is in linear ksi order and includes positive detections only, 1:nd 

 For signal-noise detectors (*detect==12) the positions nd+1 : 2nd are noise measurements 

 The integer array start (dim nc x ss x nk or nc x ss) holds the index for each isk to
   the first corresponding detection in 'signal' 

 It is possible in principle for there to be more than one detection per isk; these will
   follow in sequence, hence the name 'start' 

*/

{
    /* controls on included code (vary in pwuniform) */
    int switch0 = 1;
    int switch1 = 1;

    /* indices */
    int    i,n,g,k,m,c,s,x;
    int    wi, wxi;
    int    gi;

    /* miscellaneous */
    double temp, tempg, tempsum, tempp, templog, prwi;

    /* group arrays */
    int    *ng;       /* number per group */
    double *sumD;     
    double *sumDp;    

    /* generalised detection and detector functions */
    gfnptr gfn;
    prwfnptr prwfn;

    /* numbers of individuals and detectors */
    int    nc1;
    int    cumk[1001];   /* where this max from? */
    int    nk = 0;
    /* int    nn; */

    /* number of 'g' (detection) parameters */
    int    gpar = 2;    

    /* pre-computed detection probabilities */
    double *gk;
    double *gk0;

    /* mixture membership probability */
    double *pmix;

    /* passing data to prwfn */
    double *detspec;
    int    nval;
    int    *start = NULL;
    /* number of detections */
    int    nd = 0;  

    /* total hazard computation 2011-11-15*/
    int    *hc0 = NULL;        
    int    *hindex = NULL;     
    double *h = NULL;          
    int    next = 0;
    int    hi;
    int    fullns = 0;
    int    c0 = 0;

    /* CL only */
    int    indiv = 0;     /* indicate detection varying between individuals */
    double asum[maxnmix];
    double pd;
    int    nested = 0;  /* 0 = not nested, 1 = nested (CL only) */
    int    cumng = 0;   /* used for 'nested' */
    double dp;          /* Poisson density (cuerate) */

    double pdt = 0;

    /*-------------------------------------------------------*/

    /* MAINLINE */

     clock_t ticks1, ticks2;
     int timing = 0;
     ticks1 = clock();


/* WAS    #include "Isecrloglik.def" */
/**********************************************************************/
/*                     code shared with pwuniform                     */
/**********************************************************************/

/*
This code is included in both the 'secrloglik' and 'pwuniform' 
functions in file 'secr.c'.  'pwuniform' is used by fxi.secr.

The integer variables switch0 and switch1 control whether some 
code is executed when included in 'pwuniform':

switch0 = 0  Don't bother with 'naive' detection (always in pwuniform)
switch1 = 0  Don't bother with code used only to normalize pdf

SHARING ABANDONED BECAUSE INCLUDE FILES CONFUSE CRAN...
JUST DUPLICATE

*/

    double lambda;
    double hk;
    double par[4];   /* passing parameter values to integr fn  */
    double stdint = 1;

    /* generic failure code; reset to 0 at end                 */
    *resultcode = 1; 

    /*---------------------------------------------------------*/
    /* Adjust for zero animals                                 */
    /* (otherwise fails with nc = 0 because *nc also used to   */
    /* dimension arrays etc. so replace nc1 = max(1, *nc)      */
    if (*nc == 0) 
        nc1 = 1;
    else
        nc1 = *nc;

    /* for some purposes replaces individuals with groups     
    2011-11-21 redundant ?
    if (*like != 1) 
        nn = *gg;
    else
        nn = nc1;
    */
    /*---------------------------------------------------------*/


    /*---------------------------------------------------------*/
    ng = (int *) R_alloc(*gg, sizeof(int));

    /*---------------------------------------------------------*/
    /* Under development 2011-01-11                            */
    /* *gg > 1 with CL is used for clustered model             */
    if ((*gg > 1) && (*like == 1)) nested = 1;                 
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* Count number per group (not used for CL)                */
    /* Assume histories sorted by group = individual           */
    /* CH are numbered 0 <= n < *nc in C code                  */
    for (g=0; g<*gg; g++)
        ng[g] = 0;
    for (n=0; n<*nc; n++) { 
        g = grp[n]-1;
        ng[g]++;
    }
    /*---------------------------------------------------------*/
    /* determine number of polygons if polygon detector */
    /* for polygon detectors, kk is vector ending in zero */
    if ((*detect==3) || (*detect==4) || (*detect==6) || (*detect==7)) {
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
    }
    else
        nk = *kk;
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* dynamically allocate memory                             */
    /* use R_alloc for robust exit on interrupt                */
    /* S_alloc zeros as well 2011-11-15                        */
    /* h for total hazard added 2011-11-15                     */

    if ((*detect==3) || (*detect==4))
        nval = 2 + *cc + nc1 * *ss;
    else if ((*detect==5) || (*detect==9))
        nval = 4 + nc1 * *ss * nk;
    else if (*detect==12)    /* signalnoise */
        nval = 6 + nc1 * *ss * nk;
    else if ((*detect==6) || (*detect==7))
        nval = 2 + *cc + nc1 * *ss * nk;
    else if (*detect==8)
        nval = nc1 * *ss * nk;
    else
        nval = 4;    /* 1-4, mostly not used */

    if (timing) {
	ticks2 = clock();
	Rprintf("check 1: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

    gk = (double *) S_alloc(*cc * nk * *mm, sizeof(double));
    pmix = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    detspec = (double *) R_alloc(nval, sizeof(double));
    if (switch0)
    gk0 = (double *) S_alloc(*cc0 * nk * *mm, sizeof(double));
    if (switch1 && ((*detect==0) || (*detect==3) || (*detect==4))) {
        hc0 = (int *) R_alloc (*cc, sizeof(int));
        hindex = (int *) S_alloc (nc1 * *ss, sizeof(int));
	/*        h = (double *) S_alloc (*cc * *mm * *nmix, sizeof(double)); */
	h = (double *) S_alloc (nc1 * *ss * *mm * *nmix, sizeof(double)); 
    }
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* Identify start positions of ancillary data for each     */
    /* animal                                                  */

    /* max one detection per occasion */
    if ( (*detect == 3) || (*detect==4) ) {
        /* start[z] indexes the row in xy
           for each detection z, where z is w-order (is) */
        start = (int *) R_alloc(nc1 * *ss, sizeof(int));
        for (s=0; s< *ss; s++) {
            for (i=0; i< *nc; i++) {
                wi = *nc * s + i;
                start[wi] = nd;
                nd += (w[wi] != 0);
            }
        }
    }
    if (((*detect>=5) && (*detect<=9)) || (*detect==12)) {
        /* start[z] indexes the first row in xy (or element in signal)
           for each possible count z, where z is w-order (isk) */
        start = (int *) R_alloc(nc1 * *ss * nk, sizeof(int));
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
    /*---------------------------------------------------------*/

    /* retrieve detection function  (getgfn in utils.c)        */
    gfn = getgfn(*fn);

    /*---------------------------------------------------------*/
    /* mixture proportions                                     */
    /* by group for full likelihood, by animal otherwise       */

    for (i=0; i < nc1 * *nmix; i++) pmix[i] = 1; /* default */
    if (*nmix>1) {
        /* one extra real parameter */
        gpar++;

        for (n=0; n<nc1; n++) {
            for (x=0; x<*nmix; x++) {
                wxi = i4(n,0,0,x,*nc,*ss,nk);
                c = PIA[wxi] - 1;
                if (*like != 1) {
                    g = grp[n]-1;
                    pmix[*nmix * g + x] = gsbval[*cc * (gpar-1) + c];
                }
                else
                    pmix[*nmix * n + x] = gsbval[*cc * (gpar-1) + c];
            }
        }
    }

    /*---------------------------------------------------------*/
    /*---------------------------------------------------------*/
    /* populate pre-computed arrays                            */
    /*
        *detect may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  exclusive polygon detector
        4  exclusive transect detector
        5  signal detectors
        6  polygon detector
        7  transect detector
        8  times
        9  cue
	12 signalnoise
    */
    /*---------------------------------------------------------*/

    /* all these detectors are Bernoulli, or similar           */
    /* so we override binomN                                   */
    if ((*detect==0) || (*detect==1) || (*detect==3) || 
        (*detect==4) || (*detect==5) || (*detect==12)) {
        *binomN = 1;
    }
    if (timing) {
	ticks2 = clock();
	Rprintf("check 2: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

    if (switch1) {
    if (*detect == 0) {
        for (k=0; k<*kk; k++) {
            for (m=0; m<*mm; m++) {
                if (switch0)
                for (c=0; c<*cc0; c++) {
                    gi = i3(c,k,m,*cc0,nk);
                    gk0[gi] = gfn(k, m, c, gsb0val,
                        *cc0, traps, mask, nk, *mm, miscparm);
                }
                for (c=0; c<*cc; c++) {
                    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = gfn(k, m, c, gsbval,
                        *cc, traps, mask, nk, *mm, miscparm);
                }
            }
        }
    }
    else if ((*detect == 1) || (*detect == 2) || (*detect==5) || (*detect==8)
	     || (*detect==9) || (*detect==12)) {
        for (k=0; k<*kk; k++) {
            for (m=0; m<*mm; m++) {
                if (switch0)
                for (c=0; c<*cc0; c++) {
                    lambda = gfn(k, m, c, gsb0val, *cc0, traps, mask, *kk,
                        *mm, miscparm);
                    gi = i3(c,k,m,*cc0,nk);
                    gk0[gi] = 1 - countp(0, *binomN, lambda);
                }
                for (c=0; c<*cc; c++) {
                    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = gfn(k, m, c, gsbval, *cc,
                        traps, mask, nk, *mm, miscparm);
                }
            }
        }
    }
    else if ((*detect == 3) || (*detect == 6)) {
        if (switch0)
        for (c=0; c<*cc0; c++) {
            par[0] = gsb0val[c];
            par[1] = gsb0val[*cc0 + c];
            par[2] = gsb0val[2* *cc0 + c];
            stdint = gintegral(*fn, par);
            for (k=0; k<nk; k++) {               /* over parts */
                for (m=0; m<*mm; m++) {
                    hk = par[0] * integral2D (*fn, m, 0, par, 1,
                        traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
                    gi = i3(c, k, m, *cc0, nk);
                    gk0[gi] = 1 - countp (0, *binomN, hk);
                }
            }
        }
        R_CheckUserInterrupt();
        for (c=0; c<*cc; c++) {
            par[0] = gsbval[c];
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];
            stdint = gintegral(*fn, par);
            detspec[2+c] = stdint;               /* passed to prwipolygon */
            for (k=0; k<nk; k++) {               /* over parts */
                for (m=0; m<*mm; m++) {
                    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = par[0] * integral2D (*fn, m, 0, par, 1, traps, mask,
                        cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
                }
            }
        }
    }
    else if ((*detect == 4) || (*detect == 7)) {
        if (switch0)
        for (c=0; c<*cc0; c++) {
            par[0] = gsb0val[c];
            par[1] = gsb0val[*cc0 + c];
            par[2] = gsb0val[2* *cc0 + c];
            stdint = gintegral1(*fn, par);
            for (k=0; k<nk; k++) {               /* over transects */
                for (m=0; m<*mm; m++) {
                    hk = par[0] * integral1D (*fn, m, c, gsb0val, *cc0, traps, mask,
                        cumk[k], cumk[k+1]-1, cumk[nk], *mm)/stdint;
                    gi = i3(c,k,m,*cc0,nk);
                    gk0[gi] = 1 - countp (0, *binomN, hk);
                }
            }
        }
        R_CheckUserInterrupt();
        for (c=0; c<*cc; c++) {
            par[0] = gsbval[c];
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];
            stdint = gintegral1(*fn, par);
            detspec[2+c] = stdint;               /* passed to prwitransect */
            for (k=0; k<nk; k++) {               /* over transects */
                for (m=0; m<*mm; m++) {
                    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = par[0] * integral1D (*fn, m, c, gsbval, *cc,
                        traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
                }
            }
        }
    }
    else error ("unrecognised detector type in external C fn secrloglik");
    }  /* end switch1 */
    R_CheckUserInterrupt();

    /*=================================================================*/

    prwfn = prwicount;   /* default */
    if ((*detect == 0) && (*param == 0))
        prwfn = prwimulti;
    if ((*detect == 0) && (*param == 1))
        prwfn = prwimultiGR;
    if (*detect == 3)
        prwfn = prwipolygonX;
    if (*detect == 4)
        prwfn = prwitransectX;
    if ((*detect == 5) || (*detect==9)) 
	prwfn = prwisignal;
    if (*detect == 6)
        prwfn = prwipolygon;
    if (*detect == 7)
        prwfn = prwitransect;
    if (*detect == 8)
        prwfn = prwitimes;
    if (*detect == 12) {
	prwfn = prwisignalnoise;   /* experimental 2012-02-07 */
    }

    if (timing) {
	ticks2 = clock();
	Rprintf("check 3: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

    /*=================================================================*/
    /* need nk, pmix and gk0 for pndot() in case there are no captures */
    if (*nc == 0) goto eval;

    /*=================================================================*/
    /* detector-specific data to pass to prwi functions */

    /* total hazard for animal n on occasion s wrt mask point m */
    /* construct index 'hindex' to values in 'h'                */
    /* 'bk' model results in within-trap variation dependent on */
    /* n.s, so h must be recalc each time                       */
    /* c0 -- index of parameters for trap 0, mixture 0          */
    /* hc0[c0] -- maps c0 to sequential index 'next'            */
    /* next -- new index of parameters for each n,s             */
    /* h -- array of computed hazard for [m,next]               */

    /* mixtures are group-specific for full likelihood, and     */
    /* individual-specific for conditional likelihood           */

    if (switch1 &&((*detect == 0) || (*detect == 3) || (*detect == 4))) {

        /* recognise when not fully specified by n.s, c  */
        /* this arises when model = bk, Bk               */
        /* and when detector covariates vary by time     */
        fullns = 0;
        for (n=0; n < nc1; n++) {
            for (s=0; s < *ss; s++) {
		for (k=1; k<nk; k++) {
                   if (PIA[i4(n,s,k,0, nc1, *ss, nk)] != PIA[i4(n,s,0,0, nc1, *ss, nk)]) {
                       fullns = 1;
		       break;
		   }               
                } 
                if (fullns == 1) break;
            }
            if (fullns == 1) break;
        }

	if (timing) {
	    ticks2 = clock();
	    Rprintf("check 3.2: time used %15d ticks\n", ticks2-ticks1); 
	    ticks1 = ticks2;
	}

        for (i=0; i<*cc; i++) hc0[i] = -1;
        next = 0;        
        for (n=0; n < nc1; n++) {
            for (s=0; s < *ss; s++) {
                hi = s*nc1 + n;
                /* Case 1. within-trap variation */
		if (fullns) {
                    for (m=0; m< *mm; m++) { 
                        for (x = 0; x < *nmix; x++) {
                            for (k=0; k < nk; k++) {
                               c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
                               if (c >= 0) {
                                   gi = i3(c,k,m,*cc,nk);
                                   h[i3(x,m,hi,*nmix, *mm)] += hazard(gk[gi]);
			       }
                            }
                        }
                    }
                    hindex[hi] = hi;   
		}
                /* Case 2. no within-trap variation */
                else {
                   c0 = PIA[i4(n,s,0,0, nc1, *ss, nk)] - 1;                    
                   if (hc0[c0] < 0) {
                        hc0[c0] = next;
                        next ++;
                        for (m=0; m< *mm; m++) { 
                            for (x = 0; x < *nmix; x++) {
                                for (k=0; k < nk; k++) {
                                    c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
                                    if (c >= 0) {
                                        gi = i3(c,k,m,*cc,nk);
                                        h[i3(x,m, hc0[c0],*nmix, *mm)] += hazard(gk[gi]);
    			            }
                                }
                            }
                        }
    		    }
                    hindex[hi] = hc0[c0];
		}
            }
        }
    }
    if (timing) {
	ticks2 = clock();
	Rprintf("check 3.4: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

    if ((*detect == 3) || (*detect == 4)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        for (i=0; i< (*nc * *ss); i++)
            detspec[2+*cc+i] = (double) start[i];
    }
    else if ((*detect == 5) || (*detect==9)) {
        for (i=0; i<3; i++) detspec[i]= miscparm[i];
        detspec[3]= ((*fn == 11) || (*fn == 13));     /* spherical */
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[4+i] = (double) start[i];
    }
    else if (*detect==12) {                           /* signal-noise */        
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;                     /* number of detectors */
	detspec[2]= miscparm[0];                      /* cut */
	detspec[3]= miscparm[1];                      /* noise mean */
	detspec[4]= miscparm[2];                      /* noise sd */
        detspec[5]= ((*fn == 11) || (*fn == 13));     /* spherical */
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[6+i] = (double) start[i];
    }
    else if ((*detect == 6) || (*detect == 7)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[2+*cc+i] = (double) start[i];
    }
    else if (*detect == 8) {
        for (i=0; i< (*nc* *ss * nk); i++)
            detspec[i] = (double) start[i];
    }
    /*===============================================================*/

eval:    /* skip to here from within include block if no captures */

/* end of code shared with pwuniform */
/**********************************************************************/

    /*-------------------------------------------------------*/
    sumD = (double *) R_alloc(*gg, sizeof(double));
    sumDp = (double *) R_alloc(*gg, sizeof(double));
    /*-------------------------------------------------------*/

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
                    wxi = i4(0,s,k,x,*nc,*ss,nk);        /* changed from *kk 2011-02-07 */
                    i = PIA0[wxi];
                    for (n=1; n<*nc; n++) {
                        wxi = i4(n,s,k,x,*nc,*ss,nk);    /* changed from *kk 2011-02-07 */
                        if (i != PIA0[wxi]) {
                            indiv = 1; break;
                        }
                    }
                }
            }
        }
	if (indiv == 0) {
            /* all individuals the same */
            /* save time by doing this once, rather than inside n loop */
            for (x=0; x<*nmix; x++) {
                asum[x] = 0;
                for (m=0; m<*mm; m++) {
                    if (nested) {
                        pd = pndot (m, 0, 1, *ss, x, *nc, PIA0, gk0, *ss, nk,
                            *cc0, *nmix, gsb0val, *param);
                        asum[x] += pndotgrp (pd, miscparm[0], 0.0001);
                    }
                    else {
			    asum[x] += pndot (m, 0, 1, *ss, x, *nc, PIA0, gk0, *ss, nk,
					      *cc0, *nmix, gsb0val, *param);
                    }
                }
            }
        }
        /* else asum calculated for each individual in loop below */
        *value = 0;

        if (nested) {   /* 2011-01-11 cuerate model - specialised */

            /* Loop over individuals... */
            for (g=0; g<*gg; g++) {
                tempsum = 0;
                for (m=0; m<*mm; m++) {
                    tempg = 1;
                    /* loop over cues for this individual */
                    for (i = 0; i<ng[g]; i++) {
                        n = cumng + i;
			prwi = prwfn (m, n, 1, *ss, 0, w, xy, signal, PIA, gk,
					*binomN, detspec, h, hindex, *cc, *nc, nk, 
					*ss, *mm, 1, gfn, gsbval, traps, mask, *minprob);
			tempg *= prwi;
                    }
                    tempsum += tempg;
                }
                a[g] = asum[0];
                dp = dpois(ng[g], miscparm[0], 1);
                templog = log(tempsum) - log(a[g]) + dp;
                a[g] = *area * a[g];

                if (!R_FINITE(templog)) *resultcode = 9;
		if (*resultcode == 9) return;

                *value += templog;
                R_CheckUserInterrupt();
                cumng += ng[g];
            }     /* end loop over groups (= individuals) */
        }
        else {   /* not nested (standard option) */

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
			    asum[x] += pndot (m, n, 1, *ss, x, *nc, PIA0, gk0,
					      *ss, nk, *cc0, *nmix, gsb0val, *param);
			prwi = prwfn (m, n, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, 
				      detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix,
				      gfn, gsbval, traps, mask, *minprob);
			temp += prwi;
                    }
                    a[n] += pmix[*nmix * n + x] * asum[x];
                    tempsum += pmix[*nmix * n + x] * temp;
                }    /* end of loop over mixtures */

                templog = log(tempsum/a[n]);
                a[n] = *area * a[n];

                if (!R_FINITE(templog)) *resultcode = 9;
		if (*resultcode == 9) return;
		*value += templog;
                R_CheckUserInterrupt();
            }        /* end of loop over individuals */
        }
    }
    /*-------------------------------------------------------------------------------------------*/

    else {  /* *like==0,2  Full or Partial likelihood */

    if (timing) {
	ticks2 = clock();
	Rprintf("check 3.5: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

        for (g=0; g<*gg; g++) {
            sumD[g] = 0;
            sumDp[g] = 0;
            for (m=0; m<*mm; m++)  {
                sumD[g] += Dmask[*mm * g + m];
                for (x=0; x<*nmix; x++) {
		    pdt = pndot (m, g, 1, *ss, x, *gg, PIA0, gk0,
				 *ss, nk, *cc0, *nmix, gsb0val, *param);
                    sumDp[g] += pdt * pmix[*nmix * g + x] * Dmask[*mm * g + m];
                }
            }
        }
    if (timing) {
	ticks2 = clock();
	Rprintf("check 4: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

        *value = 0;
        /* compute likelihood component from pr(wi) */
        if (*like == 0)   /* Full likelihood only */
        for (n=0; n < *nc; n++) {
            g = grp[n]-1;
            temp = 0;
            for (x=0; x<*nmix; x++) {
                for (m=0; m<*mm; m++) {
		    prwi = prwfn (m, n, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, 
				detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, 
				gfn, gsbval, traps, mask, *minprob);
		    tempp = prwi * pmix[*nmix * g + x] * Dmask[*mm * g + m];	    
                    /* Simply aborting at this point does not work 2011-01-30 */
		    /* so following condition is not used */
                    if (tempp < 0) {
                        *resultcode = 8;
                        return;
                    }
                    else
                        temp += tempp;
                }
            }
            templog = log(temp);
            if (!R_FINITE(templog)) *resultcode = 9;
	    if (*resultcode == 9) return;
            *value += templog;

            R_CheckUserInterrupt();
        }
        /* compute likelihood component due to n */
        for (g=0; g<*gg; g++) {
            *value -= ng[g] * log(sumDp[g]);
            /* Poisson */
            if (*distrib==0) *value += gpois(ng[g], sumDp[g] * *area, 1);
            /* binomial */
            if (*distrib==1) *value += gbinomFP (ng[g], sumD[g] * *area, sumDp[g] / sumD[g], 1);
            /* superbinomial (specified N) */
            if (*distrib>=2) *value += gbinomFP (ng[g], (double) *distrib,
                sumDp[g] * *area / *distrib, 1);
        }
    }
    if (timing) {
	ticks2 = clock();
	Rprintf("check 5: time used %15d ticks\n", ticks2-ticks1); 
	ticks1 = ticks2;
    }

    *resultcode = 0;   /* successful termination secrloglik */
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

    if (*fn != 0)
    error ("invalid detection function in naived");

    for (n=0; n<*nc; n++)
    {
        x = animals[n];
        y = animals[n + *nc];

        for (i=0; i<*kk; i++)
            for (j=0; j<(i-1); j++) {

                dij = (traps[i] - traps[j]) * (traps[i] - traps[j]) +
                        (traps[i+*kk] - traps[j+*kk]) * (traps[i+*kk] - traps[j+*kk]);
                d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+*kk] - y) * (traps[i+*kk] - y);
                d22 = (traps[j] - x) * (traps[j] - x) + (traps[j+*kk] - y) * (traps[j+*kk] - y);

                if ((d21<=truncate2) && (d22<=truncate2))
                   p1p2 = exp(-(d21+d22) / 2 / *sigma / *sigma);
                else
                   p1p2 = 0;

                sump  += p1p2;
                sumdp += p1p2 * sqrt(dij);

            }
        for (i=0; i<*kk; i++) {  /* diagonal */
            d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+*kk] - y) * (traps[i+*kk] - y);
            if (d21<=truncate2)                                     /* d21=d22 */
                sump += exp(-2*d21 /2 / *sigma / *sigma)/2;
        }
    }
    *value = sumdp/sump;
}
/*==============================================================================*/

/* 2010-11-01
  replaced local eval of detection function with call to pfn
  increased truncate2 to 9 sigma^2
  added z argument
  allowed fn 0-8
*/
void naiveRPSV (
  double *sigma,   /* Parameter : detection scale */
  double *z,       /* parameter : detection shape (probably fixed) */
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
    double miscparm [3];

    if (*fn > 8)
    error ("invalid detection function in naiveRPSV");

    for (n=0; n<*nc; n++)
    {
        x = animals[n];
        y = animals[n + *nc];
        pdot   = 1;
        sumd2k = 0;
        sumpk  = 0;
        for (k=0; k< *kk; k++) {
            d2k = (traps[k] - x) * (traps[k] - x) + (traps[k + *kk]-y) * (traps[k + *kk]-y);
	    pk = pfn(*fn, d2k, 1.0, *sigma, *z, miscparm, truncate2);
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

    if (*fn != 0)
    error ("invalid detection function in naivecap2");

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

/* code for fxi individual range centre pdf */
void pwuniform (
    int    *which,       /* which one: 0 <= which < *nc (0<=which<*gg for cues)*/
    int    *xx,          /* number of points */
    double *X,           /* points at which to evaluate Pr(wi||X) */
    int    *like,        /* likelihood 0 full, 1 conditional */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
    int    *param,       /* parameterisation 0 Borchers & Efford 1 Gardner & Royle */
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

    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *PIA,         /* lookup which g0/sigma/b combination to use for given n, S, K */

    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameters (cuerate etc.) */
    int    *normal,      /* code 0 don't normalise, 1 normalise */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
/*    double *cut,       transformed signal strength threshold for detection */
    double *minprob,     /* minimum value of P(detection history) */
    double *value,       /* return values */
    int    *resultcode   /* 0 if OK */
)
{
    int    i,j,n,g,k,m,c,s,x;
    int    *ng;       /* number per group */
    int    wi;
    int    wxi;
    int    gi;
    double *pmix;
    double temp;
    double sumprwi = 1.0;
    double prwi;

    double *gk;
    double *gkx;
    double *detspec;
    double *detspecx;

    double *gk0;
    int    *hc0 = NULL;        
    int    *hindex = NULL;     
    double *h = NULL;          

    int    *hc0x;
    int    *hindexx;
    double *hx;            /* 2011-11-15 */ 
    int    *start = NULL;
    int    nval;

    gfnptr gfn;
    prwfnptr prwfn;
    int cumk[1001];
    int nk = 0;
    int nd = 0;
    int nc1 = 0;
    /* int nn; */
    int gpar = 2;    /* number of 'g' (detection) parameters */

    int nested = 0;  /* 0 = not nested, 1 = nested (CL only) */
    int *cumng;  /* used for 'nested' */
    double tempg;

    int    next = 0;
    int    hi;
    int    fullns = 0;
    int    c0 = 0;
    double p;

    int switch0;   /* compute gk0 etc.? */
    int switch1;   /* normalise or not? */

    /* dummy pointers to keep included code sweet */
    double *gsb0val;     
    int    *cc0;     

    /*===============================================================*/

    /* MAINLINE */

    /*-----------------*/
    switch0 = 0;
    switch1 = *normal;
    /*-----------------*/

/* WAS    #include "Isecrloglik.def" */
/**********************************************************************/
/*                     code shared with secrloglik                    */
/**********************************************************************/
/*

This code is included in both the 'secrloglik' and 'pwuniform' 
functions in file 'secr.c'.  'pwuniform' is used by fxi.secr.

The integer variables switch0 and switch1 control whether some 
code is executed when included in 'pwuniform':

switch0 = 0  Don't bother with 'naive' detection (always in pwuniform)
switch1 = 0  Don't bother with code used only to normalize pdf

*/

    double lambda;
    double hk;
    double par[4];   /* passing parameter values to integr fn  */
    double stdint = 1;

    /* generic failure code; reset to 0 at end                 */
    *resultcode = 1; 

    /*---------------------------------------------------------*/
    /* Adjust for zero animals                                 */
    /* (otherwise fails with nc = 0 because *nc also used to   */
    /* dimension arrays etc. so replace nc1 = max(1, *nc)      */
    if (*nc == 0) 
        nc1 = 1;
    else
        nc1 = *nc;

    /* for some purposes replaces individuals with groups     
    2011-11-21 redundant ?
    if (*like != 1) 
        nn = *gg;
    else
        nn = nc1;
    */

    /*---------------------------------------------------------*/


    /*---------------------------------------------------------*/
    ng = (int *) R_alloc(*gg, sizeof(int));

    /*---------------------------------------------------------*/
    /* Under development 2011-01-11                            */
    /* *gg > 1 with CL is used for clustered model             */
    if ((*gg > 1) && (*like == 1)) nested = 1;                 
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* Count number per group (not used for CL)                */
    /* Assume histories sorted by group = individual           */
    /* CH are numbered 0 <= n < *nc in C code                  */
    for (g=0; g<*gg; g++)
        ng[g] = 0;
    for (n=0; n<*nc; n++) { 
        g = grp[n]-1;
        ng[g]++;
    }
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* determine number of polygons if polygon detector */
    /* for polygon detectors, kk is vector ending in zero */
    if ((*detect==3) || (*detect==4) || (*detect==6) || (*detect==7)) {
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
    }
    else
        nk = *kk;
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* dynamically allocate memory                             */
    /* use R_alloc for robust exit on interrupt                */
    /* S_alloc zeros as well 2011-11-15                        */
    /* h for total hazard added 2011-11-15                     */

    if ((*detect==3) || (*detect==4))
        nval = 2 + *cc + nc1 * *ss;
    else if ((*detect==5) || (*detect==9))
        nval = 4 + nc1 * *ss * nk;
    else if ((*detect==6) || (*detect==7))
        nval = 2 + *cc + nc1 * *ss * nk;
    else if (*detect==8)
        nval = nc1 * *ss * nk;
    else
        nval = 4;    /* 1-4, mostly not used */

    gk = (double *) S_alloc(*cc * nk * *mm, sizeof(double));
    pmix = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    detspec = (double *) R_alloc(nval, sizeof(double));
    if (switch0)
    gk0 = (double *) S_alloc(*cc0 * nk * *mm, sizeof(double));
    if (switch1 && ((*detect==0) || (*detect==3) || (*detect==4))) {
        hc0 = (int *) R_alloc (*cc, sizeof(int));
        hindex = (int *) S_alloc (nc1 * *ss, sizeof(int));
	/* h = (double *) S_alloc (*cc * *mm * *nmix, sizeof(double)); */
	h = (double *) S_alloc (nc1 * *ss * *mm * *nmix, sizeof(double)); 
    }
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* Identify start positions of ancillary data for each     */
    /* animal                                                  */

    /* max one detection per occasion */
    if ( (*detect == 3) || (*detect==4) ) {
        /* start[z] indexes the row in xy
           for each detection z, where z is w-order (is) */
        start = (int *) R_alloc(nc1 * *ss, sizeof(int));
        for (s=0; s< *ss; s++) {
            for (i=0; i< *nc; i++) {
                wi = *nc * s + i;
                start[wi] = nd;
                nd += (w[wi] != 0);
            }
        }
    }
    if ((*detect>=5) && (*detect<=9)) {
        /* start[z] indexes the first row in xy (or element in signal)
           for each possible count z, where z is w-order (isk) */
        start = (int *) R_alloc(nc1 * *ss * nk, sizeof(int));
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
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* retrieve detection function  (getgfn in detectfn.c)        */
    gfn = getgfn(*fn);
    /*---------------------------------------------------------*/

    /*---------------------------------------------------------*/
    /* mixture proportions                                     */
    /* by group for full likelihood, by animal otherwise       */

    for (i=0; i < nc1 * *nmix; i++) pmix[i] = 1; /* default */
    if (*nmix>1) {
        /* one extra real parameter */
        gpar++;

        for (n=0; n<nc1; n++) {
            for (x=0; x<*nmix; x++) {
                wxi = i4(n,0,0,x,*nc,*ss,nk);
                c = PIA[wxi] - 1;
                if (*like != 1) {
                    g = grp[n]-1;
                    pmix[*nmix * g + x] = gsbval[*cc * (gpar-1) + c];
                }
                else
                    pmix[*nmix * n + x] = gsbval[*cc * (gpar-1) + c];
            }
        }
    }

    /*---------------------------------------------------------*/
    /*---------------------------------------------------------*/
    /* populate pre-computed arrays                            */
    /*
        *detect may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  exclusive polygon detector
        4  exclusive transect detector
        5  signal detectors
        6  polygon detector
        7  transect detector
        8  times
        9  cue
    */
    /*---------------------------------------------------------*/

    /* all these detectors are Bernoulli, or similar           */
    /* so we override binomN                                   */
    if ((*detect==0) || (*detect==1) || (*detect==3) || 
        (*detect==4) || (*detect==5)) {
        *binomN = 1;
    }

    if (switch1) {
	if (*detect == 0) {
	    for (k=0; k<*kk; k++) {
		for (m=0; m<*mm; m++) {
		    if (switch0)
			for (c=0; c<*cc0; c++) {
			    gi = i3(c,k,m,*cc0,nk);
			    gk0[gi] = gfn(k, m, c, gsb0val,
					  *cc0, traps, mask, nk, *mm, miscparm);
			}
		    for (c=0; c<*cc; c++) {
			gi = i3(c,k,m,*cc,nk);
			gk[gi] = gfn(k, m, c, gsbval,
			     *cc, traps, mask, nk, *mm, miscparm);
		    }
		}
	    }
	}
	else if ((*detect == 1) || (*detect == 2) || (*detect==5) || (*detect==8) || (*detect==9)) {
	    for (k=0; k<*kk; k++) {
		for (m=0; m<*mm; m++) {
		    if (switch0)
			for (c=0; c<*cc0; c++) {
			    lambda = gfn(k, m, c, gsb0val, *cc0, traps, mask, *kk,
					 *mm, miscparm);
			    gi = i3(c,k,m,*cc0,nk);
			    gk0[gi] = 1 - countp(0, *binomN, lambda);
			}
		    for (c=0; c<*cc; c++) {
			gi = i3(c,k,m,*cc,nk);
			gk[gi] = gfn(k, m, c, gsbval, *cc,
				     traps, mask, nk, *mm, miscparm);
		    }
		}
	    }
	}
	else if ((*detect == 3) || (*detect == 6)) {
	    if (switch0)
		for (c=0; c<*cc0; c++) {
		    par[0] = gsb0val[c];
		    par[1] = gsb0val[*cc0 + c];
		    par[2] = gsb0val[2* *cc0 + c];
		    stdint = gintegral(*fn, par);
		    for (k=0; k<nk; k++) {               /* over parts */
			for (m=0; m<*mm; m++) {
			    hk = par[0] * integral2D (*fn, m, 0, par, 1,
			      traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
			    gi = i3(c, k, m, *cc0, nk);
			    gk0[gi] = 1 - countp (0, *binomN, hk);
			}
		    }
		}
	    R_CheckUserInterrupt();
	    for (c=0; c<*cc; c++) {
		par[0] = gsbval[c];
		par[1] = gsbval[*cc + c];
		par[2] = gsbval[2* *cc + c];
		stdint = gintegral(*fn, par);
		detspec[2+c] = stdint;               /* passed to prwipolygon */
		for (k=0; k<nk; k++) {               /* over parts */
		    for (m=0; m<*mm; m++) {
			gi = i3(c,k,m,*cc,nk);
			gk[gi] = par[0] * integral2D (*fn, m, 0, par, 1, traps, mask,
			      cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
		    }
		}
	    }
	}
	else if ((*detect == 4) || (*detect == 7)) {
	    if (switch0)
		for (c=0; c<*cc0; c++) {
		    par[0] = gsb0val[c];
		    par[1] = gsb0val[*cc0 + c];
		    par[2] = gsb0val[2* *cc0 + c];
		    stdint = gintegral1(*fn, par);
		    for (k=0; k<nk; k++) {               /* over transects */
			for (m=0; m<*mm; m++) {
			    hk = par[0] * integral1D (*fn, m, c, gsb0val, *cc0, traps, mask,
						      cumk[k], cumk[k+1]-1, cumk[nk], *mm)/stdint;
			    gi = i3(c,k,m,*cc0,nk);
			    gk0[gi] = 1 - countp (0, *binomN, hk);
			}
		    }
		}
	    R_CheckUserInterrupt();
	    for (c=0; c<*cc; c++) {
		par[0] = gsbval[c];
		par[1] = gsbval[*cc + c];
		par[2] = gsbval[2* *cc + c];
		stdint = gintegral1(*fn, par);
		detspec[2+c] = stdint;               /* passed to prwitransect */
		for (k=0; k<nk; k++) {               /* over transects */
		    for (m=0; m<*mm; m++) {
			gi = i3(c,k,m,*cc,nk);
			gk[gi] = par[0] * integral1D (*fn, m, c, gsbval, *cc,
			      traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm) / stdint;
		    }
		}
	    }
	}
	else error ("unrecognised detector type in external C fn secrloglik");
    }  /* end switch1 */
    R_CheckUserInterrupt();

    /*=================================================================*/

    prwfn = prwicount;   /* default */
    if ((*detect == 0) && (*param == 0))
        prwfn = prwimulti;
    if ((*detect == 0) && (*param == 1))
        prwfn = prwimultiGR;
    if (*detect == 3)
        prwfn = prwipolygonX;
    if (*detect == 4)
        prwfn = prwitransectX;
    if ((*detect == 5) || (*detect==9))
        prwfn = prwisignal;
    if (*detect == 6)
        prwfn = prwipolygon;
    if (*detect == 7)
        prwfn = prwitransect;
    if (*detect == 8)
        prwfn = prwitimes;
    if (*detect == 12) {
	prwfn = prwisignalnoise;   /* experimental 2012-02-07 */
    }

    /*=================================================================*/
    /* need nk, pmix and gk0 for pndot() in case there are no captures */
    if (*nc == 0) goto eval;

    /*=================================================================*/
    /* detector-specific data to pass to prwi functions */

    /* total hazard for animal n on occasion s wrt mask point m */
    /* construct index 'hindex' to values in 'h'                */
    /* c0 -- index of parameters for trap 0, mixture 0          */
    /* hc0[c0] -- maps c0 to sequential index 'next'            */
    /* next -- new index of parameters for each n,s             */
    /* h -- array of computed hazard for [m,next]               */

    /* mixtures are group-specific for full likelihood, and     */
    /* individual-specific for conditional likelihood           */

    if (switch1 &&((*detect == 0) || (*detect == 3) || (*detect == 4))) {

        /* recognise when not fully specified by n.s, c  */
        /* this arises when model = bk, Bk               */
        fullns = 0;
        for (n=0; n < nc1; n++) {
            for (s=0; s < *ss; s++) {
		for (k=1; k<nk; k++) {
                   if (PIA[i4(n,s,k,0, nc1, *ss, nk)] != PIA[i4(n,s,0,0, nc1, *ss, nk)]) {
                       fullns = 1;
		       break;
		   }               
                } 
                if (fullns == 1) break;
            }
            if (fullns == 1) break;
        }

        for (i=0; i<*cc; i++) hc0[i] = -1;
        next = 0;        
        for (n=0; n < nc1; n++) {
            for (s=0; s < *ss; s++) {
                hi = s*nc1 + n;
                /* Case 1. within-trap variation */
		if (fullns) {
                    for (m=0; m< *mm; m++) { 
                        for (x = 0; x < *nmix; x++) {
                            for (k=0; k < nk; k++) {
                               c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
                               if (c >= 0) {
                                   gi = i3(c,k,m,*cc,nk);
                                   h[i3(x,m,hi,*nmix, *mm)] += hazard(gk[gi]);
			       }
                            }
                        }
                    }
                    hindex[hi] = hi;   
		}
                /* Case 2. no within-trap variation */
                else {
                   c0 = PIA[i4(n,s,0,0, nc1, *ss, nk)] - 1;                    
                   if (hc0[c0] < 0) {
                        hc0[c0] = next;
                        next ++;
                        for (m=0; m< *mm; m++) { 
                            for (x = 0; x < *nmix; x++) {
                                for (k=0; k < nk; k++) {
                                    c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
                                    if (c >= 0) {
                                        gi = i3(c,k,m,*cc,nk);
                                        h[i3(x,m, hc0[c0],*nmix, *mm)] += hazard(gk[gi]);
    			            }
                                }
                            }
                        }
    		    }
                    hindex[hi] = hc0[c0];
		}
            }
        }
    }

    if ((*detect == 3) || (*detect == 4)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        for (i=0; i< (*nc * *ss); i++)
            detspec[2+*cc+i] = (double) start[i];
    }
    else if ((*detect == 5) || (*detect==9)) {
        for (i=0; i<3; i++) detspec[i]= miscparm[i];
        detspec[3]= (*fn == 11);     /* spherical */
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[4+i] = (double) start[i];
    }
    else if ((*detect == 6) || (*detect == 7)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[2+*cc+i] = (double) start[i];
    }
    else if (*detect == 8) {
        for (i=0; i< (*nc* *ss * nk); i++)
            detspec[i] = (double) start[i];
    }
    /*===============================================================*/

eval:    /* skip to here from within include block if no captures */

/**********************************************************************/
/*               end of code shared with secrloglik                   */
/**********************************************************************/

    /*--------------------------------------------------------*/

    cumng = (int *) R_alloc(*gg, sizeof(int));
    cumng[0] = 0;
    for (g=1; g<*gg; g++)
       cumng[g] = cumng[g-1] + ng[g-1];

    if (nested) {
        if (*nmix>1)
            error ("cue rate models may not be combined with mixtures");
        if ((*which-1) > *gg)
            error ("requested individual exceeds number available");
    }
    else {
        if (*gg > 1)
            error("fxi does not allow for groups");
    }

    if (*normal > 0) {
	/* Rprintf("normalising\n"); */
        sumprwi = 0;
        for (x=0; x<*nmix; x++) {
            temp = 0;
            for (m=0; m<*mm; m++) {
		prwi = prwfn (m, *which-1, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, detspec,
                    h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, gfn, gsbval, traps, mask, *minprob);
		temp += prwi;
            }
            sumprwi += pmix[x] * temp;
        }    /* end of loop over mixtures */
        if (sumprwi< fuzz)
            error("zero prwi in external function pwuniform");
    }

    /*---------------------------------------------------------*/
    /* dynamically allocate memory                             */
    /* use R_alloc for robust exit on interrupt                */
    /* S_alloc zeros as well 2011-11-15                        */
    /* h for total hazard added 2011-11-15                     */

    gkx = (double *) S_alloc(*cc * nk * *xx, sizeof(double));
    detspecx = (double *) R_alloc(nval, sizeof(double));
    if ((*detect==0) || (*detect==3) || (*detect==4)) {
        hc0x = (int *) R_alloc (*cc, sizeof(int));
        hindexx = (int *) S_alloc (nc1 * *ss, sizeof(int));
        hx = (double *) S_alloc (*cc * *xx, sizeof(double));
    }
    else {
        hc0x = (int *) R_alloc (1, sizeof(int));
        hindexx = (int *) S_alloc (1, sizeof(int));
        hx = (double *) R_alloc (1, sizeof(double));
    }


    /* RECOMPUTE gk for requested X */

    if ((*detect == 0) || (*detect == 1) || (*detect == 2) || (*detect == 5) ||
             (*detect == 8) || (*detect == 9)) {
        for (k=0; k < nk; k++) {
            for (m=0; m<*xx; m++) {
                for (c=0; c<*cc; c++) {
                    gi = i3(c,k,m,*cc,nk);
                    gkx[gi] = gfn(k, m, c, gsbval, *cc,
                        traps, X, nk, *xx, miscparm);
                }
            }
        }
    }
    else if ((*detect == 3) || (*detect == 6)) {
        for (c=0; c<*cc; c++) {
            par[0] = gsbval[c];
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];
            stdint = gintegral(*fn, par);
            detspecx[2+c] = stdint;               /* passed to prwipolygon */
            for (k=0; k<nk; k++) {               /* over parts */
                for (m=0; m<*xx; m++) {
                    gi = i3(c,k,m,*cc,nk);
                    gkx[gi] = par[0] * integral2D (*fn, m, 0, par, 1, traps, X,
                        cumk[k], cumk[k+1]-1, cumk[nk], *xx) / stdint;
                }
            }
        }
    }
    else if ((*detect == 4) || (*detect == 7)) {
        for (c=0; c<*cc; c++) {
            par[0] = gsbval[c];
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];
            stdint = gintegral1(*fn, par);
            detspecx[2+c] = stdint;               /* passed to prwitransect */
            for (k=0; k<nk; k++) {               /* over transects */
                for (m=0; m<*xx; m++) {
                    gi = i3(c,k,m,*cc,nk);
                    gkx[gi] = par[0] * integral1D (*fn, m, c, gsbval, *cc,
                        traps, X, cumk[k], cumk[k+1]-1, cumk[nk], *xx) / stdint;
                }
            }
        }
    }

    /* =====================*/
    /* update detspec and h */
    /* for requested X */

    if ((*detect == 0) || (*detect == 3) || (*detect == 4)) {
        for (i=0; i<*cc; i++) hc0x[i] = -1;
        next = 0;        
        for (n=0; n < nc1; n++) {
            if (*like != 1) 
                g = grp[n]-1;
            else
                g = n;
            for (s=0; s < *ss; s++) {
               c0 = PIA[i4(n,s,0,0, nc1, *ss, nk)] - 1;
               if (hc0x[c0] < 0) {
                    hc0x[c0] = next;
                    next ++;
                    for (m=0; m< *xx; m++) { 
                        for (k=0; k < nk; k++) {
                           p = 0;
                           for (x = 0; x < *nmix; x++) {
                               c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
                               if (c >= 0) {
                                   gi = i3(c,k,m,*cc,nk);
                                   p += pmix[*nmix * g + x] * gkx[gi];
			       }
                            }
                            hx[hc0x[c0] * *xx + m] += hazard (p);
                        }
                    }
		}    
                hindexx[s*nc1 + n] = hc0x[c0];         
            }
        }
    }

    if ((*detect == 3) || (*detect == 4)) {
        detspecx[0] = (double) nk;
        detspecx[1] = (double) nd;
        for (i=0; i< (*nc* *ss); i++)
            detspecx[2+*cc+i] = (double) start[i];
    }
    else if ((*detect == 5) || (*detect==9)) {
        for (i=0; i<3; i++) detspecx[i]= miscparm[i];
        detspecx[3]= (*fn == 11) || (*fn == 12);     /* spherical */
        for (i=0; i< (*nc* *ss * nk); i++)
            detspecx[4+i] = (double) start[i];
    }
    else if ((*detect == 6) || (*detect == 7)) {
        detspecx[0] = (double) nk;
        detspecx[1] = (double) nd;
        for (i=0; i< (*nc* *ss * nk); i++)
            detspecx[2+*cc+i] = (double) start[i];
    }
    else if (*detect == 8) {
        for (i=0; i< (*nc* *ss * nk); i++)
            detspecx[i] = (double) start[i];
    }

    /* finish detspec */
    /* ===============*/

    R_CheckUserInterrupt();

    /* i indexes points at which to evaluate pdf */
    for (i=0; i< *xx; i++) {
        temp = 0;
        if (nested) {    /* no mixture */
            tempg = dpois(ng[*which-1], miscparm[0], 0);
            /* loop over cues for this individual */
            for (j = 0; j<ng[*which-1]; j++) {
                n = cumng[*which-1] + j;
		prwi = prwfn (i, n, 1, *ss, 0, w, xy, signal, PIA, gkx,
			      *binomN, detspecx, hx, hindexx, *cc, *nc, nk, *ss, *xx, 1, gfn,
			      gsbval, traps, X, *minprob);
		tempg *= prwi;
            }
            temp += tempg;
        }
        else {
            for (x=0; x<*nmix; x++) {
                temp += pmix[x] * prwfn (i, *which-1, 1, *ss, x, w, xy, signal, PIA,
/*2012-11-01  	    gkx, *binomN, detspecx, hx, hindex, *cc, *nc, nk,  */
	   	    gkx, *binomN, detspecx, hx, hindexx, *cc, *nc, nk, 
                    *ss, *xx, *nmix, gfn, gsbval, traps, X, *minprob);
            }
        }
        value[i] = temp / sumprwi;
    }

    *resultcode = 0;   /* successful termination pwuniform */
}
/*==============================================================================*/

