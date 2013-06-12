/* constants */

#define fuzz 1e-200
#define huge 1e10
#define maxnpoly 1000   
#define maxnmix 2    
#define maxvertices 200    

/* source to include */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>       /* random numbers */
#include <Rmath.h>   /* R math functions e.g. dbinom, dpois */
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>

/* This round() fn not used any more 2011-10-01 */
/* #define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5)) */

/*-------------------*/
/* data structures   */
/*-------------------*/
struct trap_animal {
    int     trap;
    int     animal;
    double  time;
};
struct rpoint {
    double x;
    double y;
};

/*-------------------*/
/* function pointers */
/*-------------------*/

typedef double (*gfnptr)(int, int, int, double[], int, double[], double[], int, int, double[]);
/*
    (int k, int m, int c, int double gsbval[], int cc, double traps[],
     double mask[], int kk, int mm, double miscparm []);
*/

typedef double (*prwfnptr)(int, int, int, int, int, int[], double[], double[], int[],
   double[], int, double[], double[], int[], int, int, int, int, int, int, gfnptr, double[],
   double[], double [], double[], double);
/*
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[], int gsb[],
     double gk[], int binomN, double detspec[], double h[], double hindex[], int cc, 
     int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[], 
     double traps[], double Tsk[], double mask[], double minp);
*/

typedef double (*fnptr)(double[], double);

/*--------------------------*/
/* functions to be exported */
/*--------------------------*/

/*---------*/
/* secr.c */
/*---------*/

void pdotpoint (double *xy, int *nxy, double *traps, int *detect, double *Tsk, int *kk, 
               int *fn, double *par, int *nocc, double *w2, int *binomN, double *value);

void pdotpoly (double *xy, int *nxy, double *traps, int *detect, double *Tsk, int *nk, 
	       int *kk, int *fn, double *par, int *nocc, int *binomN, double *value);

/*
void pdottransect (double *xy, int *nxy, double *traps, int *detect, double *Tsk, int *nk, 
               int *kk, int *fn, double *par, int *nocc, int *binomN, double *value);
*/

void integralprw1 (
    int    *detect,    /* detector 0 multi, 1 proximity etc. */
    int    *param,     /* parameterisation 0 Borchers & Efford 1 Gardner & Royle */
    double *gsb0val,   /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *grp,         /* group number for 0<=n<*nc   [full likelihood only] !!!!!*/
    int    *nc,        /* number of individuals */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    int    *gg,          /* number of groups   !!!!! */
    int    *nmix,      /* number of mixtures */
    int    *knownclass,  /* known membership of 'latent' classes */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *Tsk,       /* s x k matrix effort on occasion s at at detector k */
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
    );

/*---------------------------------------------------------------------*/

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
    int    *knownclass,  /* known membership of 'latent' classes */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *Tsk,         /* s x k matrix effort on occasion s at at detector k */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *Dmask,       /* density at each point on mask, possibly x group */
    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    double *gsb0val,     /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *cc0,         /* number of g0/sigma/b combinations for naive animals */
    int    *PIA,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *PIA0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive] */
    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameters */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
/*    double *cut,       transformed signal strength threshold for detection */
    double *minprob,     /* minimum value of P(detection history) */
    double *a,           /* a(theta) */
    double *value,       /* return value integral */
    int    *resultcode   /* 0 if OK */
    );
/*---------------------------------------------------------------------*/

void MRsecrloglik (
    int    *detect,      /* detector 0 multi, 1 proximity */
    int    *distrib,     /* distribution 0 Poisson, 1 binomial */
    int    *w,           /* capture histories (1:nc, 1:s, 1:k) */
    int    *Tu,          /* unmarked counts s = qq+1 ... ss, dim=c(*kk, *ss, *gg) */
    int    *Tm,          /* marked not ID counts s = qq+1 ... ss, dim=c(*kk, *ss, *gg) */
    int    *grp,         /* group number for 0<=n<*nc   [full likelihood only] */
    int    *nc,          /* number of capture histories */
    int    *ss,          /* total number of occasions */
    int    *qq,          /* number of marking occasions */
    int    *kk,          /* number of traps */
    int    *mm,          /* number of points on mask */
    int    *gg,          /* number of groups */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *Dmask,       /* density at each point on mask, possibly x group */
    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    double *gsb0val,     /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *cc0,         /* number of g0/sigma/b combinations for naive animals */
    int    *gsb,         /* lookup which g0/sigma/b to use for given n, S, K */
    int    *gsb0,        /* lookup which g0/sigma/b to use for given g, S, K [naive animal] */
    double *pID,         /* Parameter value - probability marked animal identified on resighting */
    double *area,        /* area associated with each mask point (ha) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard rate, 2 = exponential etc. */
    int    *binomN,      /*  */
    double *minprob,     /* minimum value of P(detection history) */
    double *value,       /* return value integral of pr(w0) */
    int    *resultcode   /* 0 if OK */
);
/*---------------------------------------------------------------------*/

void naived (
    double *sigma,      /* Parameter : detection scale */
    int    *kk,         /* number of traps */
    int    *nc,
    int    *wt,         /* integer weights */
    double *traps,      /* x,y locations of traps (first x, then y) */
    double *animals,    /* x,y locations of traps (first x, then y) */
    int    *fn,         /* code 0 = halfnormal ONLY */
    double *value       /* return value  */
);
/*---------------------------------------------------------------------*/

void naiveRPSV (
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* parameter : detection shape (probably fixed) */
    int    *kk,         /* number of traps */
    int    *nc,
    int    *wt,         /* integer weights */
    double *traps,      /* x,y locations of traps (first x, then y)   */
    double *animals,    /* x,y locations of animals (first x, then y) */
    int    *fn,         /* code 0 = halfnormal ONLY */
    double *value       /* return value */
);
/*---------------------------------------------------------------------*/

void naivecap2 (
    int    *detect,     /* code 0 = multicatch, 1 = proximity */
    double *g0,         /* Parameter : detection magnitude */
    double *sigma,      /* Parameter : detection scale */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *mm,
    int    *wt,         /* integer trap weights */
    double *traps,      /* x,y locations of traps (first x, then y) */
    double *mask,       /* x,y points on mask (first x, then y) */
    int    *fn,         /* code 0 = halfnormal ONLY */
    double *value       /* return value  */
);
/*---------------------------------------------------------------------*/

void makelookup (   
    double *x,          /* input matrix */
    int    *nrow,       /* input */
    int    *ncol,       /* input */
    int    *unique,     /* output number of unique rows */
    double *y,          /* output matrix of unique rows (byrow=T) */
    int    *index,      /* output lookup rows of x in y */
    int    *resultcode  /* zero if OK */
);
/*---------------------------------------------------------------------*/

void integral2Dtest
    (int *fn, int *m, int *c, double *gsbval, int *cc, double *traps, 
    double *mask, int *n1, int *n2, int *kk, int *mm, double *result);
/*---------------------------------------------------------------------*/

void pwuniform (
    int    *which,       /* which one: 0 <= which < *nc */
    int    *xx,          /* number of points */
    double *X,           /* points at which to evaluate Pr(wi|X) */
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
    int    *knownclass,  /* known membership of 'latent' classes */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *Tsk,         /* s x k matrix effort on occasion s at at detector k */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *gsb,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameter */
    int    *normal,      /* code 0 don't normalise, 1 normalise */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
/*    double *cut,       transformed signal strength threshold for detection */
    double *minprob,     /* minimum value of P(detection history) */
    double *value,       /* return values */
    int    *resultcode   /* 0 if OK */
    );

/*-----------*/
/* simsecr.c */
/*-----------*/

void simsecr (
    int    *detect,     /* detector 0 multi, 1 proximity */
    double *gsb0val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    double *gsb1val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [caught before] */
    int    *cc0,        /* number of g0/sigma/b combinations for naive animals */
    int    *cc1,        /* number of g0/sigma/b combinations for caught before */
    int    *gsb0,       /* lookup which g0/sigma/b to use for given g, S, K [naive animal] */
    int    *gsb1,       /* lookup which g0/sigma/b to use for given n, S, K [caught before] */
    int    *N,          /* number of animals */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *nmix,       /* number of classes */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    double *Tsk,        /* ss x kk array of 0/1 usage codes or effort */
    int    *btype,      /* code for behavioural response 0 none 1 b 2 bk 3 k */
    int    *Markov,     /* code 0 if behavioural response is learned, 1 if Markov */
    int    *binomN,     /* number of trials for 'count' detector modelled with binomial */
    double *miscparm,   /* detection threshold on transformed scale, etc. */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 = uniform */
    int    *maxperpoly, /* */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    double *detectedXY, /* x,y locations of detections  */
    double *signal,     /* vector of signal strengths, one per detection */
    int    *value,      /* return value array of trap locations n x s */
    int    *resultcode
);

/*-------------*/
/* trapping.c */
/*-------------*/

void trappingpolygon (
    double *lambda,      /* Parameter : expected detection events per hectare */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *npoly,       /* number of different polygons */
    int    *kk,          /* number of vertices + 1 (assuming closed) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    double *Tsk,         /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *binomN,      /* 0 poisson, 1 Bernoulli, or number of trials for 
                            'count' detector modelled with binomial */
    int    *maxperpoly,
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/
void trappingtelemetry (
    double *lambda,  /* Parameter : expected detection events per hectare */
    double *sigma,   /* Parameter : detection scale */
    double *z,       /* Parameter : detection shape (hazard) */
    int    *ss,      /* number of occasions */
    int    *N,       /* number of animals */
    double *animals, /* x,y points of animal range centres (first x, then y)  */
    int    *fn,      /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,      /* truncation radius */
    int    *binomN,  /* 0 poisson, 1 Bernoulli, or number of binomial trials */
    int    *maxperpoly, /*   */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
    );
/*---------------------------------------------------------------------*/

void trappingpolygonX (
    double *g0,          /* Parameter : detection intercept */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *npoly,       /* number of different polygons */
    int    *kk,          /* number of vertices + 1 (assuming closed) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    double *Tsk,         /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/

void trappingtransect (
    double *lambda,      /* Parameter : expected detection events per metre */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *ntransect,   /* number of different transects */
    int    *kk,          /* number of vertices + 1 (assuming closed) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    double *Tsk,         /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *binomN,      /* 0 poisson, 1 Bernoulli, or number of trials for 
                            'count' detector modelled with binomial */
    int    *maxperpoly,
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/

void trappingtransectX (
    double *g0,          /* Parameter : detecttion intercept */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *ntransect,   /* number of different transects */
    int    *kk,          /* number of vertices + 1 (assuming closed) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    double *Tsk,         /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/

void trappingsignal (
    double *beta0,     /* Parameter : intercept */
    double *beta1,     /* Parameter : slope */
    double *sdS,       /* Parameter : error sd */
    double *cut,       /* detection threshold on transformed scale, etc. */
    double *muN,       /* noise mean */
    double *sdN,       /* noise sd */
    double *sdM,       /* movement between occasions */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    double *Tsk,       /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,        /* code 10 = signal strength, 11 = signal strength with spherical spr */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *signal,    /* signal strength, one per detection */  
    double *noise,     /* noise, one per detection, if signalnoise */  
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/

void trappingtimes (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    double *Tsk,       /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,        /* truncation radius */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *times,     /* time of detection within occasion */  
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/

void trappingcount (
    double *g0,         /* Parameter : detection intercept */
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    double *Tsk,        /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 11 = normal signal */
    double *w2,         /* truncation radius */
    int    *binomN,     /* number of trials for 'count' detector modelled with binomial */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode
);
/*---------------------------------------------------------------------*/

void trappingmulti (
    double *g0,         /* Parameter : detection magnitude */
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    double *Tsk,        /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,         /* truncation radius */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode  /* 0 for successful completion */
);
/*---------------------------------------------------------------------*/

void trappingsingle (
    double *g0,         /* Parameter : detection magnitude */
    double *sigma,      /* Parameter : detection scale */
    double *b,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    double *Tsk,        /* ss x npoly array of 0/1 usage codes or effort */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,         /* truncation radius */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode  /* 0 for successful completion */
);
/*---------------------------------------------------------------------*/

/*------------*/
/* unmarked.c */
/*------------*/

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
    );

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
    int    *type,        /* code 0 = simple, 1 = integrated */
    double *value,       /* return value */
    int    *resultcode   /* 0 if OK */
    );

/*------------*/
/* detectfn.c */
/*------------*/

gfnptr getgfn (int fn); 
fnptr gethfn (int fn); 

/* detection functions used internally in C code */
/* these exist in two forms 'g' and 'h' */

double hn (double param [], double r);
double hr (double param [], double r);
double he (double param [], double r);
double hnc (double param [], double r);
double un (double param [], double r);
double hf (double param [], double r);
double hann (double param [], double r);
double hcln (double param [], double r);
double hcg (double param [], double r);
double hsigbin  (double param [], double r);
double hsig (double param [], double r);
double hsigsph (double param [], double r);
double hhn (double param [], double r);
double hhr (double param [], double r);
double hex (double param [], double r);
double han (double param [], double r);
double hcumg (double param [], double r);

double ghn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double ghnc 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double ghr 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double ghe 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gun
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double ghf 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gan 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gcln
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gcn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gcg
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double grs
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gsigbin
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gsig
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gsigsph
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gsigSN
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double gsigsphSN
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double lhn
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double lhr
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double lex 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double lan 
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);
double lcg
    (int k, int m, int c, double gsbval[], int cc, double traps[],
    double mask[], int kk, int mm, double miscparm[]);

double pfn (
    int fn,
    double d2val,
    double g0,
    double sigma,
    double z,
    double miscparm[],
    double w2);

/*---------*/
/* utils.c */
/*---------*/

int i3 (
    int i, 
    int j, 
    int k, 
    int ii, 
    int jj);
/*---------------------------------------------------------------------*/

int i4 (
    int i, 
    int j, 
    int k, 
    int l, 
    int ii, 
    int jj, 
    int kk);
/*---------------------------------------------------------------------*/

double d2 (
    int k,
    int m,
    double A1[],
    double A2[],
    int A1rows,
    int A2rows);

/*--------------------------------------------------------------------------*/

/* customised dpois */
double gpois (int count, double lambda, int uselog);
/*--------------------------------------------------------------------------*/

/* customised dbinom */
double gbinom(int count, int size, double p, int uselog);
/*--------------------------------------------------------------------------*/

/* customised dnbinom parameterised as size, mu */
double gnbinom (int count, int size, double mu, int uselog);
/*--------------------------------------------------------------------------*/

/* binomial density allowing non-integer (floating point) size */
double gbinomFP (int count, double size, double p, int uselog);
/*--------------------------------------------------------------------------*/

/* probability of count with distribution specified by binomN */
double countp (int count, int binomN, double lambda);
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
    int spherical);
/*---------------------------------------------------------------------*/

void SegCircle (
    double *p1x, double *p1y, 
    double *p2x, double *p2y, 
    double *scx, double *scy, 
    double *r, 
    double *seg);
/*---------------------------------------------------------------------*/

double SegCircle2 (
    double p1x, double p1y, 
    double p2x, double p2y, 
    double scx, double scy, 
    double r);
/*---------------------------------------------------------------------*/

double expmin (double x);

double Random ();

double distance (struct rpoint p1, struct rpoint p2);

double rcount (int binomN, double lambda, double Tsk);
/*---------------------------------------------------------------------*/

struct rpoint getxy(
    double l, 
    double cumd[], 
    struct rpoint line[], 
    int kk, 
    double offset);
/*---------------------------------------------------------------------*/

double gintegral (int fn, double par[]);
/*---------------------------------------------------------------------*/

double gintegral1 (
    int fn, 
    double par[]);
/*---------------------------------------------------------------------*/

double randomtime (double p);
/*---------------------------------------------------------------------*/

void probsort (
    int n, 
    struct trap_animal tran[]);
/*---------------------------------------------------------------------*/

double gr (
    int *fn, 
    double par[], 
    struct rpoint xy, 
    struct rpoint animal);
/*---------------------------------------------------------------------*/

/* random point from 2-D radial distribution specified by g function */
void gxy (
    int *n, 
    int *fn,
    double *par,
    double *w, 
    double *xy
); 
/*---------------------------------------------------------------------*/

void nearest (   
    int    *nxy,        /* number of input points */
    double *xy,         /* input points */
    int    *ntrap,      /* input */
    double *traps,      /* input */
    int    *p,          /* output index of nearest point */
    double *d           /* output distance to nearest point */
);
/*---------------------------------------------------------------------*/

void inside (
    double *xy, 
    int *n1, 
    int *n2, 
    int *np, 
    double *poly, 
    int *in);
/*---------------------------------------------------------------------*/

void ontransect (
    double *xy, 
    int    *n1, 
    int    *n2, 
    int    *npts, 
    double *transect,
    double *tol, 
    int    *on);
/*---------------------------------------------------------------------*/

void alongtransect (
    double *xy,
    int    *n1,
    int    *n2,
    int    *npts,
    double *transect,
    double *tol,
    double *along);
/*---------------------------------------------------------------------*/

/* for transect detectors */
double integral1D
    (int fn, int m, int c, double gsbval[], int cc, double traps[],
     double mask[], int n1, int n2, int kk, int mm, double ex[]);
/*---------------------------------------------------------------------*/

/* for polygon detectors */
double integral2D  (int fn, int m, int c, double gsbval[], int cc, double traps[],
		    double mask[], int n1, int n2, int kk, int mm, double ex[]);

/*---------------------------------------------------------------------*/

