#include "secr.h"

double minimumexp = -100;

/*--------------------------------------------------------------------------*/

double expmin (double x)
{
  if (x < minimumexp)
      return(0);
  else
      return(exp(x));
}
/*--------------------------------------------------------------------------*/

/* index to vector element corresponding to cell i,j,k in 3D array
   stored in column-major order */

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}
/*--------------------------------------------------------------------------*/

/* index to vector element corresponding to cell i,j,k,l in 4D array
   stored in column-major order */

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}
/*--------------------------------------------------------------------------*/

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
/*--------------------------------------------------------------------------*/
/* customised dpois */
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
/*--------------------------------------------------------------------------*/

/* customised dbinom */
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
/*--------------------------------------------------------------------------*/

/* customised dnbinom parameterised as size, mu */
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
/*--------------------------------------------------------------------------*/

/* binomial density allowing non-integer (floating point) size */
double gbinomFP (int count, double size, double p, int uselog)
{
    return ( lgamma(size+1) - lgamma(size-count+1) - lgamma(count+1) +
             count * log(p) + (size - count) * log (1-p) );
}
/*--------------------------------------------------------------------------*/

/* probability of count with distribution specified by binomN */
double countp (int count, int binomN, double lambda) {
    /* Poisson */
    if (binomN == 0) {
	if (count == 0) 
            return (exp(-lambda));
	else
	    return (dpois(count, lambda, 0));
        /* return ( gpois (count, lambda, 0)); replaced 2012-12-18 */
    }

    /* Bernoulli */
    else if (binomN == 1) {
        if (count == 0)
            return ( 1 - lambda );
        else
            return ( lambda );
    }

    /* negative binomial */
    else if (binomN < 0)
        return ( gnbinom (count, binomN, lambda, 0) );

    /* binomial */
    else
        return ( gbinom (count, binomN, lambda, 0) );
/*        return ( gbinom (count, binomN, lambda / binomN, 0) ); replaced 2012-12-23 */
}
/*--------------------------------------------------------------------------*/

double distance (struct rpoint p1, struct rpoint p2) {
    return(sqrt ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y)));
}
/*--------------------------------------------------------------------------*/

/* random point from 2-D radial distribution specified by g function */
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
/*--------------------------------------------------------------------------*/

void nearest (
  int    *nxy,      /* number of input points */
  double *xy,       /* input points */
  int    *ntrap,    /* input */
  double *traps,    /* input */
  int    *p,        /* output indices of nearest point */
  double *d)        /* output distances to nearest point */
{
    int i,j;
    int id=-1;
    double d2;
    double d2min;
    for (j=0; j<*nxy; j++) {
        id = -1;
        d2min = 1e100;
        for (i=0; i<*ntrap; i++)
        {
            d2 = (traps[i] - xy[j]) * (traps[i] - xy[j]) +
                 (traps[i + *ntrap] - xy[j+ *nxy]) * (traps[i + *ntrap] - xy[j + *nxy]);
            if (d2 < d2min) { d2min = d2; id = i; }
        }
        d[j] = sqrt(d2min);
        p[j] = id+1;
    }
}
/*--------------------------------------------------------------------------*/

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
        if (fabs(d)>0) { N = N/fabs(d);  d = d/fabs(d); }
        theta += atan2(N, d);
    }
    theta = fabs(theta);
    if (fabs(theta - 2* M_PI) < cutoff)    /* M_PI is Rmath.h constant */
        *in = 1;    /* inside */
    else
        *in = 0;    /* outside */
}
/*--------------------------------------------------------------------------*/

struct rpoint getxy(double l, double cumd[], struct rpoint line[], 
    int kk, double offset) {
/* return the xy coordinates of point l metres along a transect */
/* offset is the starting position for this transect */
    int k;
    double pr, d, d12;
    struct rpoint xy;
    for (k=offset+1; k<(offset+kk); k++) 
        if (cumd[k]>l) break;
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
/*--------------------------------------------------------------------------*/

double mufn (
    int k,
    int m,
    double b0,
    double b1,
    double A1[],
    double A2[],
    int A1rows,
    int A2rows,
    int spherical)
/*
   Return predicted signal strength at m for source at point k,
   given strength at source of b0 dB and attenuation of b1 dB/m.
   Spherical spreading is included if spherical > 0
   Coordinates of points are in A1 and A2 which have respectively
   A1rows and A2rows
*/
{
    double d2val;
    d2val = d2(k,m, A1, A2, A1rows, A2rows);
    if (spherical <= 0)
	return (b0 + b1 * sqrt(d2val));
    else {
	if (d2val>1) {
	    return (b0 - 10 * log ( d2val ) / 2.302585 + b1 * (sqrt(d2val)-1)); 
	}
	else
	    return (b0);
    }

}
/*==============================================================================*/

/* Calculate the length of intersection of a line segment and a circle
   Based on C code of Paul Bourke November 1992
   Line segment is defined from p1 to p2
   Circle is of radius r and centred at sc
   Two potential points of intersection given by
     p = p1 + mu1 (p2-p1)
     p = p1 + mu2 (p2-p1)
   Return 0 if line segment does not intersect circle
*/

void SegCircle (
    double *p1x, double *p1y, 
    double *p2x, double *p2y, 
    double *scx, double *scy, 
    double *r, 
    double *seg) 
{
    double a,b,c;
    double bb4ac;
    double dpx;
    double dpy;

    double mu1;
    double mu2;
    int p1in;
    int p2in;
 
    double i1x;
    double i1y;
    double i2x;
    double i2y;

    int i1between;
    double d1,d2;

    *seg = 0;

    /* case where both p1 and p2 inside circle */

    /*
      Rprintf ("p1 %6.3f %6.3f\n", *p1x, *p1y);
      Rprintf ("p2 %6.3f %6.3f\n", *p2x, *p2y);
      Rprintf ("sc %6.3f %6.3f\n", *scx, *scy);
      Rprintf ("r %6.3f \n", *r);
    */

    p1in = ((*scx - *p1x) * (*scx - *p1x) + 
	    (*scy - *p1y) * (*scy - *p1y)) < (*r * *r);
    p2in = ((*scx - *p2x) * (*scx - *p2x) + 
	    (*scy - *p2y) * (*scy - *p2y)) < (*r * *r);
    dpx = *p2x - *p1x;
    dpy = *p2y - *p1y;

    a = dpx * dpx + dpy * dpy;
    b = 2 * (dpx * (*p1x - *scx) + dpy * (*p1y - *scy));
    c = *scx * *scx + *scy * *scy;
    c += *p1x * *p1x + *p1y * *p1y;
    c -= 2 * (*scx * *p1x + *scy * *p1y);
    c -= *r * *r;
    bb4ac = b * b - 4 * a * c;

    /* case of no intersection */
    if ((fabs(a) < 1e-10) || (bb4ac < 0)) {
        *seg = 0;   
    }
    else {

        mu1 = (-b + sqrt(bb4ac)) / (2 * a);
        mu2 = (-b - sqrt(bb4ac)) / (2 * a);

        i1x = *p1x + mu1 * (*p2x - *p1x);
        i1y = *p1y + mu1 * (*p2y - *p1y);
        i2x = *p1x + mu2 * (*p2x - *p1x);
        i2y = *p1y + mu2 * (*p2y - *p1y);
        if (((mu1<0) && (mu2<0)) || ((mu1>1) && (mu2>1))) {
            /* no intersection */
            *seg = 0;
        }
        else {
            if (((mu1<0) && (mu2>1)) || ((mu1>1) && (mu2<0))) {
                /* both inside */
                *seg = sqrt ((*p1x - *p2x) * (*p1x - *p2x) + 
                    (*p1y - *p2y) * (*p1y - *p2y));
            }
            else {
                if ((mu1>0) && (mu1<1) && (mu2>0) && (mu2<1)) {
        	        /* two intersections */
                    *seg = sqrt ((i1x - i2x) * (i1x - i2x) + 
                          (i1y - i2y) * (i1y - i2y));
    	        }
                else {
                    /* one intersection */
                    d1 = sqrt((i1x - *p1x) * (i1x * *p1x) + 
                         (i1y - *p1y) * (i1y - *p1y));
                    d2 = sqrt((i1x - *p2x) * (i1x * *p2x) + 
                         (i1y - *p2y) * (i1y - *p2y));
                    i1between = sqrt(a) < (d1 + d2 + 1e-10);
                    if (p1in) {
                        if (i1between) {
                            i2x = *p1x;
                            i2y = *p1y;
                        }
                        else {
                            i1x = *p1x;
                            i1y = *p1y;
                        }
                    }
                    if (p2in) {
                        if (i1between) {
                            i2x = *p2x;
                            i2y = *p2y;
                        }
                        else {
                            i1x = *p2x;
                            i1y = *p2y;
                        }
                    }
                    *seg = sqrt ((i1x - i2x) * (i1x - i2x) + 
                        (i1y - i2y) * (i1y - i2y));
		}
	    }
        }    
    }
}

double SegCircle2 (
    double p1x, double p1y, 
    double p2x, double p2y, 
    double scx, double scy, 
    double r
    ) 
{
    double a,b,c;
    double bb4ac;
    double dpx;
    double dpy;

    double mu1;
    double mu2;
    int p1in;
    int p2in;
 
    double i1x;
    double i1y;
    double i2x;
    double i2y;

    int i1between;
    double d1,d2;
    double seg = 0;

    /* case where both p1 and p2 inside circle */

    /*
      Rprintf ("p1 %6.3f %6.3f\n", p1x, p1y);
      Rprintf ("p2 %6.3f %6.3f\n", p2x, p2y);
      Rprintf ("sc %6.3f %6.3f\n", scx, scy);
      Rprintf ("r %6.3f \n", r);
    */

    p1in = ((scx - p1x) * (scx - p1x) + 
	    (scy - p1y) * (scy - p1y)) < (r * r);
    p2in = ((scx - p2x) * (scx - p2x) + 
	    (scy - p2y) * (scy - p2y)) < (r * r);
    if (p1in && p2in) {        
        seg = sqrt ((p1x - p2x) * (p1x - p2x) + 
                     (p1y - p2y) * (p1y - p2y));
        return (seg);
    }

    dpx = p2x - p1x;
    dpy = p2y - p1y;

    a = dpx * dpx + dpy * dpy;
    b = 2 * (dpx * (p1x - scx) + dpy * (p1y - scy));
    c = scx * scx + scy * scy;
    c += p1x * p1x + p1y * p1y;
    c -= 2 * (scx * p1x + scy * p1y);
    c -= r * r;
    bb4ac = b * b - 4 * a * c;

    /* case of no intersection */
    if ((fabs(a) < 1e-10) || (bb4ac < 0)) {
        return (0);   
    }

    mu1 = (-b + sqrt(bb4ac)) / (2 * a);
    mu2 = (-b - sqrt(bb4ac)) / (2 * a);

    i1x = p1x + mu1 * (p2x - p1x);
    i1y = p1y + mu1 * (p2y - p1y);
    i2x = p1x + mu2 * (p2x - p1x);
    i2y = p1y + mu2 * (p2y - p1y);

    if (((mu1<0) && (mu2<0)) || ((mu1>1) && (mu2>1))) {
        /* no intersection */
        seg = 0;
    }
    else {
        if (((mu1<0) && (mu2>1)) || ((mu1>1) && (mu2<0))) {
            /* both inside */
	    seg = sqrt ((p1x - p2x) * (p1x - p2x) + 
                (p1y - p2y) * (p1y - p2y));
        }
        else {
            if ((mu1>0) && (mu1<1) && (mu2>0) && (mu2<1)) {
    	        /* two intersections */
                seg = sqrt ((i1x - i2x) * (i1x - i2x) + 
                      (i1y - i2y) * (i1y - i2y));
	    }
            else {
                /* one intersection */
                d1 = sqrt((i1x - p1x) * (i1x * p1x) + 
                     (i1y - p1y) * (i1y - p1y));
                d2 = sqrt((i1x - p2x) * (i1x * p2x) + 
                     (i1y - p2y) * (i1y - p2y));
                i1between = sqrt(a) < (d1 + d2 + 1e-10);
                if (p1in) {
                    if (i1between) {
                        i2x = p1x;
                        i2y = p1y;
                    }
                    else {
                        i1x = p1x;
                        i1y = p1y;
                    }
                }
                if (p2in) {
                    if (i1between) {
                        i2x = p2x;
                        i2y = p2y;
                    }
                    else {
                        i1x = p2x;
                        i1y = p2y;
                    }
                }
                seg = sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
	    }
	}
    }
    return(seg);    
}

/*----------------------------------------------------------------*/

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
/*----------------------------------------------------------------*/

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
/*----------------------------------------------------------------*/

/*
   This call to R random number generator requires preceding
       GetRNGstate();
   and following
       PutRNGstate();
*/
double Random () {
    return( unif_rand() );
}
/*----------------------------------------------------------------*/

/* random count from different distributions */

double rcount (int binomN, double lambda, double Tsk) {

    /* Poisson */
    if (binomN == 0)
        return ( rpois(lambda * Tsk) );

    /* negative binomial */
    else if (binomN < 0) {
        /* must use 'size, prob' parameters */
        /* prob = size / (size + mu) */
        binomN = abs(binomN);
        return ( rnbinom(binomN, binomN / (binomN+ (lambda * Tsk))) );
    }

    else { 
        if (fabs(Tsk-1) > 1e-10)               /* not 1.0 */
	    lambda = 1 - pow(1-lambda, Tsk);   /* 2012-12-18 */

	/* Bernoulli */
	if (binomN == 1) {
	    if (Random() < lambda)
		return (1);
	    else
		return (0);
	}

	/* binomial */
	else
/*	    return ( rbinom(binomN, lambda / binomN) ); changed 2012-12-23 */
	    return ( rbinom(binomN, lambda) );
    }
}
/*----------------------------------------------------------------*/

void rgr(double *x, int n, void *ex) {
    int i;
    int fn;
    double * p;
    double tmp[4];
    p = (double*) ex;
    for (i=0; i<4; i++) tmp[i] = p[i];
    fn = tmp[3];
    fnptr fnp = hn;
    fnp = gethfn(fn);
    for (i=0; i<n; i++) {
        x[i] = x[i] * fnp(tmp,x[i]);   /* r.g(r) */
    }
}

void justgr(double *x, int n, void *ex) {
    int i;
    int fn;
    double * p;
    double tmp[4];
    p = (double*) ex;
    for (i=0; i<4; i++) tmp[i] = p[i];
    fn = tmp[3];
    fnptr fnp = hn;
    fnp = gethfn(fn);
    for (i=0; i<n; i++) {
        x[i] = fnp(tmp,x[i]);   /* g(r) */
    }
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
/*----------------------------------------------------------------*/

double gintegral1 (int fn, double par[]) {
/* integral of 1-D function */
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

    /* treat uniform separately */
    if (fn == 4) {
        return (2 * par[0] * par[1]);  /* symmetric about 0 */
    }

    a = 0;
    b = 1;    /* signals bounds 0,Inf */
    ex[0] = par[0];
    ex[1] = par[1];
    ex[2] = par[2];
    ex[3] = fn;

    Rdqagi(justgr, &ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, 
        &ier, &limit, &lenw, &last, iwork, work);

    /* ignoring ier etc. */
    return (result * 2);    /* symmetric about 0 */
}
/*----------------------------------------------------------------*/

double gr (int *fn, double par[], struct rpoint xy, struct rpoint animal) {
    double r;
    fnptr fnp = hn;
    fnp = gethfn(*fn);
    r = distance (xy, animal);
    return (fnp(par,r));
}
/*----------------------------------------------------------------*/

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
/*----------------------------------------------------------------*/

void alongtransect (
    double *xy,
    int    *n1,
    int    *n2,
    int    *npts,
    double *transect,
    double *tol,
    double *along)
{
/*
    How far is point xy from start of transect? 2011-06-07
    We assume transect coordinates are in col-major order (x's then y's)
    http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/ 29/11/09
*/
    int k;
    double r;
    double u;
    struct rpoint p,p1,p2,p3;

    p3.x = xy[0];
    p3.y = xy[1];

    *along = 0;

    for (k= *n1; k < *n2; k++)
    {
        p1.x = transect[k];
        p1.y = transect[k+*npts];
        r = distance (p1,p3);
        if (r < *tol) {
            return;
        } 

        p2.x = transect[k+1];
        p2.y = transect[k+1+*npts];
        if (distance(p1,p2) > 0) {
            u = ((p3.x-p1.x) * (p2.x-p1.x) + (p3.y-p1.y) * (p2.y-p1.y)) /
                ((p2.x-p1.x) * (p2.x-p1.x) + (p2.y-p1.y) * (p2.y-p1.y));
            if ((u>=0) && (u<=1)) {
                p.x = p1.x + u * (p2.x-p1.x);
                p.y = p1.y + u * (p2.y-p1.y);
                r = distance (p,p3);
                if (r < *tol) {
		    *along += distance(p,p1);
                    return;
                } 
            }
            *along += distance(p1,p2);
        }
    }
}
/*----------------------------------------------------------------*/
