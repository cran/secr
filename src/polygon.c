/* C code for R package 'secr' */
/* Murray Efford */

/*

This file contains functions for integration of home range overlap
with polygon detectors and similar. Moved from secr.c 2012-11-13 

The choice of detection function is a bit problematic. Before Jan 2016 (2.10.1) 
the function integrated was g(r) (a probability, hence between 0 and 1) but the theory 
is for the hazard. Testing zfn 2016-01-02.

Shifting entirely to hazard formulation (zfn) 2017-03-22

*/

#include "secr.h"

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
    fnptr fnzr = zhnr;
    p = (double*) ex;
    fn = round(p[3]);
    mx = p[4];
    my = p[5];
    xy[0] = p[6];

    /* set detection function */
    fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
    for (i=0; i<n; i++) {
        xy[1] = x[i];   /* set each y value */
        d = sqrt ( (xy[1]-my)*(xy[1]-my) + (xy[0]-mx)*(xy[0]-mx) );
        x[i] = fnzr(p, d);   /* z(r) */
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
        yab(x, &i, &kk, poly, &a, &b);   /* refine limits here */
        p[6] = x[i];                   /* pass forward the value of x; consider &ex etc. */
        Rdqags(fy, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
          &limit, &lenw, &last, iwork, work);
        x[i] = result;
    }
}

double integral2D  (int fn, int m, int c, double gsbval[], int cc, double traps[],
		    double mask[], int n1, int n2, int kk, int mm, double ex[]) {
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

void fx1 (double *x, int n, void *ex) {
    int i;
    int ns;
    int fn;
/*    struct rpoint *line;*/
    struct rpoint line[maxvertices * 2];
    struct rpoint mxy;
    struct rpoint xy;
    double * p;
    double cumd[maxvertices * 2];
    double d;
    fnptr fnzr = zhnr;
    /* extract parameters passed in void pointer ex */
    p = (double*) ex;
    fn = round(p[3]);
    mxy.x = p[4];
    mxy.y = p[5];
    ns = round(p[9]);
    /* coordinates of vertices */
/*    line = (struct rpoint *) R_alloc(ns, sizeof(struct rpoint));*/
    for (i=0; i<ns; i++) {
        line[i].x = p[i+10];
        line[i].y = p[i+ns+10];
    }
    /* cumulative distance along line */
    /* cumd = (double *) R_alloc(ns + 1, sizeof(double)); */
    cumd[0] = 0;
    for (i=0; i<(ns-1); i++) {
        cumd[i+1] = cumd[i] + distance (line[i],line[i+1]);
    }
    /* set detection function - default zhnr */
    fnzr = getzfnr(fn);  // 2016-01-02, 2017-03-22
    /* for each x in x[] */
    for (i=0; i<n; i++) {
        xy = getxy (x[i], cumd, line, ns, 0);
        d = distance (xy, mxy);
        x[i] = fnzr(p, d);   /* z(r) */
    }
}

double integral1D
    (int fn, int m, int c, double gsbval[], int cc, double traps[],
     double mask[], int n1, int n2, int kk, int mm, double ex[])
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
