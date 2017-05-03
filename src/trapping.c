/*
    trappingXXXX routines perform simulated sampling of 2D popn with various
    detector types

trappingsingle
trappingmulti
trappingproximity
trappingcount
trappingpolygon
trappingpolygonX
trappingtransect
trappingtransectX
trappingsignal
trappingtelemetry
trappingtimes

2013-04-04 extended to allow k-specific g0
2015-12-21 bug fix in trapping count for Poisson

*/

#include "secr.h"

void trappingsingle (
    double *g0,        /* Parameter : detection magnitude  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    double *dist2,     /* distances squared (optional: -1 if unused) */
    double *Tsk,       /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 uniform */
    double *w2,        /* truncation radius */
    int    *binomN,    /* not used */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode /* 0 for successful completion */
)
{
    int    i,j,k,s,gi,si;
    int    nc         = 0;
    int    tr_an_indx = 0;
    double d2val;
    double p;
    int nanimals;         /* temporary */
    int ntraps;           /* temporary */
    int occupied[*kk];    /* today */
    int intrap[*N];       /* today   */

    struct  trap_animal *tran;
    double event_time;
    int anum = 0;
    int tnum = 0;
    int nextcombo;
    int finished;
    int OK;
    double miscparm[3];
    double Tski;

    /* MAIN LINE */
    *resultcode = 1;
    GetRNGstate();
    tran = (struct trap_animal *) R_alloc(*N * *kk, sizeof(struct trap_animal));

    /* ------------------------------------------------------ */
    /* pre-compute distances */
    if (dist2[0] < 0) {
	dist2 = (double *) S_alloc(*kk * *N, sizeof(double));
	makedist2 (*kk, *N, traps, animals, dist2);
    }
    else {
	squaredist (*kk, *N, dist2);
    }

    /* ------------------------------------------------------ */

    for (i=0; i<*N; i++) caught[i] = 0;   /* has animal i been caught in session? */
    for (s=0; s<*ss; s++) {
        /* initialise day */
        tr_an_indx = 0;
        nanimals = *N;
        ntraps   = *kk;
        for (i=0; i<*N; i++) intrap[i] = 0;
        for (k=0; k<*kk; k++) occupied[k] = 0;
        nextcombo = 0;

        /* make tran */
        for (i=0; i<*N; i++) {                     /* animals */
	    for (k=0; k<*kk; k++) {                /* traps */
		/* if (used[s * *kk + k]) { */
		Tski = Tsk[s * *kk + k];
		if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */
		    /* d2val = d2(i,k, animals, traps, *N, *kk); */
		    d2val = d2L(k, i, dist2, *kk);
		    gi = (caught[i]>0) * *ss * *kk + k * *ss + s;
                    si = (caught[i]>0) * *ss + s;
		    p = pfn(*fn, d2val, g0[gi], sigma[si], z[s], miscparm, *w2); 

		    if (fabs(Tski-1) > 1e-10)           /* 2012 12 26 */
			p = 1 - exp(Tski * log(1 - p));

		    event_time = randomtime(p);
		    if (event_time <= 1) {
			tran[tr_an_indx].time   = event_time;
			tran[tr_an_indx].animal = i;    /* 0..*N-1 */
			tran[tr_an_indx].trap   = k;    /* 0..*kk-1 */
			tr_an_indx++;
		    }
		}
	    }
	}

        if (tr_an_indx>0) probsort (tr_an_indx, tran);

        /* make captures */
        while ((nextcombo < tr_an_indx) && (nanimals>0) && (ntraps>0)) {
            finished = 0;
            OK       = 0;
            while ((1-finished)*(1-OK) > 0) {    /* until finished or OK */
                if (nextcombo >= (tr_an_indx)) finished = 1;  /* no more to process */
                else {
                    anum = tran[nextcombo].animal;
                    tnum = tran[nextcombo].trap;
                    OK = (1-occupied[tnum]) * (1-intrap[anum]); /* not occupied and not intrap */
                    nextcombo++;
                }
            }
            if (finished==0) {                   /* Record this capture */
                  occupied[tnum] = 1;
                  intrap[anum]   = tnum+1;       /* trap = k+1 */
                  nanimals--;
                  ntraps--;
            }
        }
        for (i=0; i<*N; i++)
        if (intrap[i]>0) {
            if (caught[i]==0) {                  /* first capture of this animal */
               nc++;
               caught[i] = nc;                   /* nc-th animal to be captured */
               for (j=0; j<*ss; j++)
                   value[*ss * (nc-1) + j] = 0;
             }
             value[*ss * (caught[i]-1) + s] = intrap[i];  /* trap = k+1 */
        }
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();

}
/*==============================================================================*/

void trappingmulti (
    double *g0,         /* Parameter : detection magnitude  */
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    double *dist2,     /* distances squared (optional: -1 if unused) */
    double *Tsk,        /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,         /* truncation radius */
    int    *binomN,     /* not used */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode  /* 0 for successful completion */
)

{
    double *h;
    double hsum[*N];
    double cump[*kk+1];
    double runif;
    int    i,j,k,s,gi,si;
    int    nc;
    double d2val;
    double p;
    double miscparm[3];
    double Tski;

    *resultcode = 1;
    cump[0] = 0;
    nc = 0;
    GetRNGstate();
/*    h = (double *) R_alloc(*N * *kk, sizeof(double)); */
    h = (double *) S_alloc(*N * *kk, sizeof(double));   /* initialise to zero */

    /* ------------------------------------------------------ */
    /* pre-compute distances */
    if (dist2[0] < 0) {
	dist2 = (double *) S_alloc(*kk * *N, sizeof(double));
	makedist2 (*kk, *N, traps, animals, dist2);
    }
    else {
	squaredist (*kk, *N, dist2);
    }
    /* ------------------------------------------------------ */

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            hsum[i] = 0;
            for (k=0; k<*kk; k++)
            {
		/* d2val = d2(i,k, animals, traps, *N, *kk); */
		d2val = d2L(k, i, dist2, *kk);
		gi = (caught[i]>0) * *ss * *kk + k * *ss + s;
                si = (caught[i]>0) * *ss + s;
		p = pfn(*fn, d2val, g0[gi], sigma[si], z[s], miscparm, *w2); 
		/* p = p * used[s * *kk + k];  zero if not used 2009 11 09 */
		Tski = Tsk[s * *kk + k];
		if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */
		    h[k * *N + i] = -Tski * log(1 - p);
		    hsum[i] += h[k * *N + i];
		}
            }

            for (k=0; k<*kk; k++) {
                cump[k+1] = cump[k] + h[k * *N + i]/hsum[i];
            }

            if (Random() < (1-exp(-hsum[i])))
            {
		if (caught[i]==0)           /* first capture of this animal */
		{
		    nc++;
		    caught[i] = nc;
		    for (j=0; j<*ss; j++)
			value[*ss * (nc-1) + j] = 0;
		}
		runif = Random();
		k = 0;
		while ((runif > cump[k]) && (k<*kk)) k++;  /* pick a trap */
		value[*ss * (caught[i]-1) + s] = k;
            }
        }
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();

}
/*==============================================================================*/

void trappingproximity (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    double *dist2,     /* distances squared (optional: -1 if unused) */
    double *Tsk,       /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,        /* truncation radius */
    int    *binomN,    /* 0 poisson, 1 Bernoulli, or number of trials for 'count'
                          detector modelled with binomial */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)
{
    double d2val;
    double theta;
    int    i,j,k,l,s;
    int    nc;
    int    count;
    double miscparm[3];
    double Tski;
    *resultcode = 1;
    nc = 0;
    GetRNGstate();

    /* ------------------------------------------------------ */
    /* pre-compute distances */
    if (dist2[0] < 0) {
	dist2 = (double *) S_alloc(*kk * *N, sizeof(double));
	makedist2 (*kk, *N, traps, animals, dist2);
    }
    else {
	squaredist (*kk, *N, dist2);
    }

    /* ------------------------------------------------------ */

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            for (k=0; k<*kk; k++) {
                /* if (used[s * *kk + k]) { */
		Tski = Tsk[s * *kk + k];
                if (fabs(Tski) > 1e-10) {          /* nonzero 2012 12 18 */
		    /* d2val = d2(i,k, animals, traps, *N, *kk); */
		    d2val = d2L(k, i, dist2, *kk);
                    /* theta = pfn(*fn, d2val, g0[s], sigma[s], z[s], miscparm, *w2); */
                    theta = pfn(*fn, d2val, g0[k * *ss + s], sigma[s], z[s], miscparm, *w2); 

                    if (theta>0) {
                        count = rcount (1, theta, Tski);
                        if (count>0)
                        {
                             /* first capture of this animal */
                             if (caught[i]==0)                  
                             {
                                 nc++;
                                 caught[i] = nc;
                                 for (j=0; j<*ss; j++)
                                   for (l=0; l<*kk; l++)
                                     value[*ss * ((nc-1) * *kk + l) + j] = 0;
                             }
                             value[*ss * ((caught[i]-1) * *kk + k) + s] = count;
                        }
                    }
                }
            }
        }
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*===========================================================================*/

void trappingcount (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    double *dist2,     /* distances squared (optional: -1 if unused) */
    double *Tsk,       /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,        /* truncation radius */
    int    *binomN,    /* 0 poisson, 1 Bernoulli, or number of trials for 'count'
                          detector modelled with binomial */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)

{
    double d2val;
    double theta;
    int    i,j,k,l,s;
    int    nc;
    int    count;
    double miscparm[3];
    double Tski;

    *resultcode = 1;
    nc = 0;
    GetRNGstate();

    /* ------------------------------------------------------ */
    /* pre-compute distances */
    if (dist2[0] < 0) {
	dist2 = (double *) S_alloc(*kk * *N, sizeof(double));
	makedist2 (*kk, *N, traps, animals, dist2);
    }
    else {
	squaredist (*kk, *N, dist2);
    }

    /* ------------------------------------------------------ */

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            for (k=0; k<*kk; k++) {
                /* if (used[s * *kk + k]) { */
		Tski = Tsk[s * *kk + k];
                if (fabs(Tski) > 1e-10) {          /* nonzero 2012 12 18 */
		    /* d2val = d2(i,k, animals, traps, *N, *kk); */
		    d2val = d2L(k, i, dist2, *kk);
                    /* theta = pfn(*fn, d2val, g0[s], sigma[s], z[s], miscparm, *w2); */
                    theta = pfn(*fn, d2val, g0[k * *ss + s], sigma[s], z[s], miscparm, *w2); 
                    if (theta>0) {
			if (binomN[s] == 1) {
			    count = rcount (round(Tski), theta, 1);
			}
			else {
                if (binomN[s] == 0)   
                    theta = hazard(theta);
                count = rcount (binomN[s], theta, Tski);
			}
                        if (count>0)
                        {
                             /* first capture of this animal */
                             if (caught[i]==0)                  
                             {
                                 nc++;
                                 caught[i] = nc;
                                 for (j=0; j<*ss; j++)
                                   for (l=0; l<*kk; l++)
                                     value[*ss * ((nc-1) * *kk + l) + j] = 0;
                             }
                             value[*ss * ((caught[i]-1) * *kk + k) + s] = count;
                        }
                    }
                }
            }
        }
    }
    *n = nc;

    *resultcode = 0;
    PutRNGstate();
}
/*===========================================================================*/

void trappingpolygon (
    double *lambda,  /* Parameter : expected detection events per hectare */
    double *sigma,   /* Parameter : detection scale */
    double *z,       /* Parameter : detection shape (hazard) */
    int    *ss,      /* number of occasions */
    int    *npoly,   /* number of different polygons */
    int    *kk,      /* number of vertices + 1 (assumed closed) each polygon */
    int    *N,       /* number of animals */
    double *animals, /* x,y points of animal range centres (first x, then y)  */
    double *traps,   /* x,y polygon vertices (first x, then y)  */
    double *Tsk,     /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,      /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,      /* truncation radius */
    int    *binomN,  /* 0 poisson, 1 Bernoulli, or number of binomial trials */
    int    *maxperpoly, /*   */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)

{
    int    i,j,k,l,s,t;
    int    nc = 0;
    int    np = 1;       /* number of points each call of gxy */
    int    nd = 0;
    int    count;
    double par[3];
    double w;
    double ws;
    int    maxdet;
    int    g=0;
    int    *gotcha;
    double xy[2];
    gotcha = &g;
    int cumk[maxnpoly+1];   /* limit maxnpoly polygons */
    int sumk;
    int n1,n2;
    double *workXY;
    int    *sortorder;
    double *sortkey;

    int J;
    double dx,dy,d;
    double Tski = 1.0;

    *resultcode = 1;
    GetRNGstate();
    if (*npoly>maxnpoly) {
        *resultcode = 2;
        return;
    }
    cumk[0] = 0;
    for (k =0; k<*npoly; k++) cumk[k+1] = cumk[k] + kk[k];
    sumk = cumk[*npoly];
    maxdet = *N * *ss * *npoly * *maxperpoly;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));

    /* find maximum distance between animal and detector vertex */
    w = 0;
    J = cumk[*npoly];
    for (i = 0; i< *N; i++) {
        for (j = 0; j < J; j++) {
            dx = animals[i] - traps[j];
            dy = animals[*N + i] - traps[J + j];
	    d = sqrt(dx*dx + dy*dy);
            if (d > w) w = d;
        }
    } 
    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        if (w > (10 * sigma[s])) 
            ws = 10 * sigma[s];
        else 
            ws = w;

        par[0] = lambda[s];
        par[1] = sigma[s];
        par[2] = z[s];
        if (lambda[s]>0) {
            for (i=0; i<*N; i++) {
                count = rcount (binomN[s], lambda[s], Tski);   /* NOT YET DEFINED !!!!*/
                /* require maximum at r=0 */
                if (*fn == 6) error ("annular normal not allowed in trappingpolygon");
                par[0] = 1;
                for (j=0; j<count; j++) {
                    gxy (&np, fn, par, &ws, xy);            /* simulate location */
                    xy[0] = xy[0] + animals[i];
                    xy[1] = xy[1] + animals[*N + i];
                    for (k=0; k<*npoly; k++) {             /* each polygon */
                        /* if (used[s * *npoly + k]) { */
			Tski = Tsk[s * *npoly + k];
			if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */
                            n1 = cumk[k];
                            n2 = cumk[k+1]-1;
                            inside(xy, &n1, &n2, &sumk, traps, gotcha);  /* assume closed */
                            if (*gotcha > 0) {
                                if (caught[i]==0) {            /* first capture of this animal */
                                    nc++;
                                    caught[i] = nc;
                                    for (t=0; t<*ss; t++)
                                        for (l=0; l<*npoly; l++)
                                            value[*ss * ((nc-1) * *npoly + l) + t] = 0;
                                }
                                nd++;
                                if (nd >= maxdet) {
                                    *resultcode = 2;           /* error */
                                    return;
                                }
                                value[*ss * ((caught[i]-1) * *npoly + k) + s]++;
                                workXY[(nd-1)*2] = xy[0];
                                workXY[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                            }
                        }
                    }
                }
            }
        }
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingtransect (
    double *lambda,      /* Parameter : expected detection events per metre? */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *ntransect,   /* number of different transects */
    int    *kk,          /* number of vertices per transect (vector) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, 
                            then y)  */
    double *traps,       /* x,y transect vertices (first x, then y)  */
    double *Tsk,         /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential 
                            etc. */
    double *w2,          /* truncation radius */
    int    *binomN,      /* 0 poisson, 1 Bernoulli, or number of trials for 
                            binomial */
    int    *maxperpoly,  /* */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)
{
    int    i,j,k,l,s,t;
    int    nc = 0;
    int    nd = 0;
    int    sumk;
    int    n1,n2;
    int    count;
    int    maxdet;
    int    gotcha;

    int    *cumk;
    double *cumd;
    struct rpoint *line;
    struct rpoint xy;
    struct rpoint animal;
    double par[3];
    double lx;
    double maxg = 0;
    double lambdak;
    double grx;
    double *workXY;
    int    *sortorder;
    double *sortkey;
    double *ex;
    double H;
    double Tski;

    *resultcode = 1;
    GetRNGstate();

    cumk = (int *) R_alloc(*ntransect+1, sizeof(int));
    cumk[0] = 0;
    for (k =0; k<*ntransect; k++)
        cumk[k+1] = cumk[k] + kk[k];
    sumk = cumk[*ntransect];
    line = (struct rpoint *) R_alloc(sumk, sizeof(struct rpoint));
    cumd = (double *) R_alloc(sumk, sizeof(double));
    maxdet = *N * *ss * *ntransect * *maxperpoly;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));

    /* coordinates of vertices */
    for (i=0; i<sumk; i++) {
        line[i].x = traps[i];
        line[i].y = traps[i+sumk];
    }

    /* cumulative distance along line */
    for (k=0; k<*ntransect; k++) {
        cumd[cumk[k]] = 0;
        for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
            cumd[i+1] = cumd[i] + distance(line[i], line[i+1]);
        }
    }

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {                            /* each occasion */
        if (lambda[s]>0) {
            for (i=0; i<*N; i++) {                     /* each animal */
                animal.x = animals[i];
                animal.y = animals[i + *N];
                for (k=0; k<*ntransect; k++) {         /* each transect */
                    /* if (used[s * *ntransect + k]) {  */
		    Tski = Tsk[s * *ntransect + k];
		    if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */
			
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        par[0] = lambda[s];
                        par[1] = sigma[s];
                        par[2] = z[s];
                        H = hintegral1(*fn, par) / par[0];

/* flaw in following: integral1D can exceed diameter */
                        lambdak = par[0] * integral1D (*fn, i, 0, par, 1, 
			     traps, animals, n1, n2, sumk, *N, ex) / H;
                        count = rcount(binomN[s], lambdak, Tski);
                        maxg = 0;
                        if (count>0) {    /* find maximum - approximate */
                            for (l=0; l<=100; l++) {
                                lx = cumd[n2] * l/100;
                                xy = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xy, animal);
/* Rprintf("xy %12.6f, %12.6f grx %12.6f\n", xy.x, xy.y, grx); */
                                maxg = fmax2(maxg, grx);
                            }
                            for (l=n1; l<=n2; l++) {
                                xy = line[l];
                                grx = gr (fn, par, xy, animal);
                                maxg = fmax2(maxg, grx);
                            }
                            maxg= 1.2 * maxg;               /* safety margin */
                            if (maxg<=0) maxg=0.0001;         /* not found */
                        }
/* error("stopping"); */
                        for (j=0; j<count; j++) {
                            gotcha = 0;
                            l = 0;
                            while (gotcha == 0) {
                                /* n2-n1 2011-02-08 here and elsewhere */
                                /* location along transect */
                                lx = Random() * (cumd[n2]-cumd[n1]);   
                                xy = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xy, animal);

                                /* rejection sampling */
                                if (*fn == 4) {
                                    if (grx > 1e-10)
                                        gotcha = 1;
				}
                                else {
                                    if (Random() < (grx/maxg))    
                                        gotcha = 1;
				}
                                l++;
                                if (l % 10000 == 0)
                                    R_CheckUserInterrupt();
                                /* give up and accept anything!!!! */
/* why why why? */
                                if (l>1e6) {
				    Rprintf ("trials exceeded 1e6 in trappingtransect\n"); 
                                    gotcha = 1;        
                                }
                            }
                            /* first capture of this animal */
                            if (caught[i]==0) {               
                                nc++;
                                caught[i] = nc;
                                for (t=0; t<*ss; t++)
                                    for (l=0; l<*ntransect; l++)
                                        value[*ss * ((nc-1) * *ntransect + l) 
                                           + t] = 0;
                            }
                            nd++;
                            if (nd >= maxdet) {
                                *resultcode = 2;
                                return;  /* error */
                            }
                            value[*ss * ((caught[i]-1) * *ntransect + k) + s]++;
                            workXY[(nd-1)*2] = xy.x;
                            workXY[(nd-1)*2+1] = xy.y;
                            sortkey[nd-1] = (double) (k * *N * *ss + s * *N +
                                caught[i]);
                        }
                    }
                }
            }
        }
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingpolygonX (
    double *g0,          /* Parameter : detection pr intercept */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *npoly,       /* number of different polygons */
    int    *kk,          /* number of vertices + 1 (assumed closed) for each polygon */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    double *Tsk,         /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)

{
    int    i,j,k,s,t;
    int    nc = 0;
    int    np = 1;       /* number of points each call of gxy */
    int    nd = 0;
    double par[3];
    double w;
    double ws;
    int    maxdet;
    int    g=0;
    int    *gotcha;
    double xy[2];
    gotcha = &g;
    int cumk[maxnpoly+1];   /* limit maxnpoly polygons */
    int sumk;
    int n1,n2;
    double *workXY;
    int    *sortorder;
    double *sortkey;

    int J;
    double dx,dy,d;
    int maybecaught;
    double Tski;

    *resultcode = 1;
    GetRNGstate();

    cumk[0] = 0;
    for (k =0; k<*npoly; k++) cumk[k+1] = cumk[k] + kk[k];
    sumk = cumk[*npoly];
    maxdet = *N * *ss;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));

    /* find maximum distance between animal and detector vertex */
    w = 0;
    J = cumk[*npoly];
    for (i = 0; i< *N; i++) {
        for (j = 0; j < J; j++) {
            dx = animals[i] - traps[j];
            dy = animals[*N + i] - traps[J + j];
	    d = sqrt(dx*dx + dy*dy);
            if (d > w) w = d;
        }
    } 
    for (i=0; i<*N; i++) caught[i] = 0;

    for (s=0; s<*ss; s++) {
        if (w > (10 * sigma[s])) 
            ws = 10 * sigma[s];
        else 
            ws = w;

        par[0] = g0[s];
        par[1] = sigma[s];
        par[2] = z[s];

        if (g0[s]>0) {
            for (i=0; i<*N; i++) {
                maybecaught = Random() < g0[s];
                /* require maximum at r=0 */
                if (*fn == 6) error ("annular normal not allowed in trappingpolygonX");

                par[0] = 1;

                /* simulate location */
                if (maybecaught) {
                    gxy (&np, fn, par, &ws, xy);
                    xy[0] = xy[0] + animals[i];
                    xy[1] = xy[1] + animals[*N + i];
                    /* which polygon, if any? */
                    for (k=0; k < *npoly; k++) {             /* each polygon */
                        /* if (used[s * *npoly + k]) { */
			Tski = Tsk[s * *npoly + k];
			if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */

                            n1 = cumk[k];
                            n2 = cumk[k+1]-1;
                            inside(xy, &n1, &n2, &sumk, traps, gotcha);  /* assume closed */
                            if (*gotcha > 0) {
                                if (caught[i]==0) {            /* first capture of this animal */
                                    nc++;
                                    caught[i] = nc;
                                    for (t=0; t<*ss; t++)
                                        value[*ss * (nc-1) + t] = 0;
                                }
                                nd++;
                                value[*ss * (caught[i]-1) + s] = k+1;
                                workXY[(nd-1)*2] = xy[0];
                                workXY[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (s * *N + caught[i]);
                                break;   /* no need to look at more poly */
                            }
			}
                    }
                }
            }
        }
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingtransectX (
    double *g0,          /* Parameter : expected detection events per metre? */
    double *sigma,       /* Parameter : detection scale */
    double *z,           /* Parameter : detection shape (hazard) */
    int    *ss,          /* number of occasions */
    int    *ntransect,   /* number of different transects */
    int    *kk,          /* number of vertices per transect (vector) */
    int    *N,           /* number of animals */
    double *animals,     /* x,y points of animal range centres (first x, then y)  */
    double *traps,       /* x,y polygon vertices (first x, then y)  */
    double *Tsk,         /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential etc. */
    double *w2,          /* truncation radius */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)
{
    int    i,k,l,s,t;
    int    nc = 0;
    int    nd = 0;
    int    sumk;
    int    n1,n2;
    int    count = 0;
    int    maxdet;
    int    gotcha;

    int    *cumk;
    double *cumd;
    struct rpoint *line;
    struct rpoint xy;
    struct rpoint animal;
    double par[3];
    double lx;
    double maxg = 0;
    double lambdak;
    double grx;
    double *workXY;
    int    *sortorder;
    double *sortkey;
    double *ex;
    double H;
    double sumhaz;
    double pks;
    double Tski;

    *resultcode = 1;
    GetRNGstate();

    cumk = (int *) R_alloc(*ntransect+1, sizeof(int));
    cumk[0] = 0;
    for (k =0; k<*ntransect; k++)
        cumk[k+1] = cumk[k] + kk[k];
    sumk = cumk[*ntransect];
    line = (struct rpoint *) R_alloc(sumk, sizeof(struct rpoint));
    cumd = (double *) R_alloc(sumk, sizeof(double));
 
    maxdet = *N * *ss;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));

    /* coordinates of vertices */
    for (i=0; i<sumk; i++) {
        line[i].x = traps[i];
        line[i].y = traps[i+sumk];
    }

    /* cumulative distance along line */
    for (k=0; k<*ntransect; k++) {
        cumd[cumk[k]] = 0;
        for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
            cumd[i+1] = cumd[i] + distance(line[i], line[i+1]);
        }
    }

/*    error ("implementation of simulations from exclusive transects has not been completed");*/

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {                               /* each occasion */
        if (g0[s]>0) {

            par[0] = g0[s];
            par[1] = sigma[s];
            par[2] = z[s];
            H = hintegral1(*fn, par);

            for (i=0; i<*N; i++) {                        /* each animal */
                animal.x = animals[i];
                animal.y = animals[i + *N];
                sumhaz = 0;
                for (k=0; k<*ntransect; k++) {            /* sum hazard */
                    /* if (used[s * *ntransect + k]) { */
		    Tski = Tsk[s * *ntransect + k];
		    if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */

                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        sumhaz += -log(1 - par[0] * integral1D (*fn, i, 0, par, 1, traps, 
			   animals, n1, n2, sumk, *N, ex) / H);
                    }
                }
                 
                for (k=0; k<*ntransect; k++) {            /* each transect */
                    /* if (used[s * *ntransect + k]) { */
		    Tski = Tsk[s * *ntransect + k];
		    if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        lambdak = par[0] * integral1D (*fn, i, 0, par, 1, traps, 
						       animals, n1, n2, sumk, *N, ex) / H;
    		        pks = (1 - exp(-sumhaz)) * (-log(1-lambdak)) / sumhaz;
                        count = Random() < pks;
                        maxg = 0;
    
                        if (count>0) {                        /* find maximum - approximate */
                            for (l=0; l<=100; l++) {
                                lx = (cumd[n2] - cumd[n1]) * l/100;    
                                xy = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xy, animal);
                                maxg = fmax2(maxg, grx);
                            }
                            for (l=n1; l<=n2; l++) {
                                xy = line[l];
                                grx = gr (fn, par, xy, animal);
                                maxg = fmax2(maxg, grx);
                            }
                            maxg= 1.2 * maxg;                 /* safety margin */
                            if (maxg<=0) maxg=0.0001;         /* not found */
    
                            gotcha = 0;
                            l = 0;
                            while (gotcha == 0) {
                                lx = Random() * (cumd[n2]-cumd[n1]);  /* location along transect */
                                xy = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xy, animal);
                                if (Random() < (grx/maxg))    /* rejection sampling */
                                    gotcha = 1;
                                l++;
                                if (l % 10000 == 0)
                                    R_CheckUserInterrupt();
                                if (l>1e8) gotcha = 1;        /* give up and accept anything!!!! */
                            }
                            if (caught[i]==0) {               /* first capture of this animal */
                                nc++;
                                caught[i] = nc;
                                for (t=0; t<*ss; t++)
                                    value[*ss * (nc-1) + t] = 0;
                            }
                            nd++;
                            value[*ss * (caught[i]-1) + s] = k+1;
                            workXY[(nd-1)*2] = xy.x;
                            workXY[(nd-1)*2+1] = xy.y;
                            sortkey[nd-1] = (double) (s * *N + caught[i]);
                        }
                        if (count>0) break;   /* no need to look further */
                    }
		}
            }
        }
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

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
    double *dist2,     /* distances squared (optional: -1 if unused) */
    double *Tsk,       /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,        /* code 10 = signal strength, 11 = signal strength with sph spread */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *signal,    /* signal strength, one per detection */
    double *noise,     /* noise, one per detection, if signalnoise */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)

/* returned signal strength (*fn==10) is on transformed scale */
/* limited to Bernoulli count model binomN = 1 */
{
    double muS;
    double signalvalue;
    double noisevalue;
    int    i,j,k,l,s;
    int    nc = 0;
    int    nd = 0;
    int    maxdet;
    double *worksignal;
    double *worknoise;
    int    *sortorder;
    double *sortkey;
    double animalss [*N*2];
    double Tski;

    *resultcode = 1;
    GetRNGstate();

    maxdet = *N * *ss * *kk;
    worksignal = (double*) R_alloc(maxdet, sizeof(double));
    worknoise = (double*) R_alloc(maxdet, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    for (i=0; i<*N; i++) {
        animalss[i] = animals[i];
        animalss[*N+i] = animals[*N+i];
    }

    /* ------------------------------------------------------ */
    /* pre-compute distances */
    if (dist2[0] < 0) {
	dist2 = (double *) S_alloc(*kk * *N, sizeof(double));
	makedist2 (*kk, *N, traps, animalss, dist2);
    }
    else {
	squaredist (*kk, *N, dist2);
    }

    /* ------------------------------------------------------ */

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            if (*sdM > fuzz) {
                animalss[i] = animals[i] +  norm_rand() * *sdM;
                animalss[*N+i] = animals[*N+i] + norm_rand() * *sdM;
            }
            else {
            }
            for (k=0; k<*kk; k++) {
                /* if (used[s * *kk + k]) { */
		Tski = Tsk[s * *kk + k];
		if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */

                    /*
                    if ((*fn == 10) || (*fn == 12))
                        muS  = mufn (i, k, beta0[s], beta1[s], animalss, traps, *N, *kk, 0);
                    else
                        muS  = mufn (i, k, beta0[s], beta1[s], animalss, traps, *N, *kk, 1);
                    */

                    if ((*fn == 10) || (*fn == 12))
                        muS  = mufnL (k, i, beta0[s], beta1[s], dist2, *kk, 0);
                    else
                        muS  = mufnL (k, i, beta0[s], beta1[s], dist2, *kk, 1);
                    signalvalue = norm_rand() * sdS[s] + muS;
                    if ((*fn == 12) || (*fn == 13)) {
			noisevalue = norm_rand() * sdN[s] + muN[s];  
			/* 2012-09-19 shouldn't this be SminusN? */
			if ((signalvalue-noisevalue) > *cut) {
			    if (caught[i]==0) {              /* first capture of this animal */
				nc++;
				caught[i] = nc;
				for (j=0; j<*ss; j++)
				    for (l=0; l<*kk; l++)
					value[*ss * ((nc-1) * *kk + l) + j] = 0;
			    }
			    nd++;
			    if (nd > maxdet) {
				*resultcode = 2;
				return;  /* error */
			    }
			    value[*ss * ((caught[i]-1) * *kk + k) + s] = 1;
			    worksignal[nd-1] = signalvalue;
			    worknoise[nd-1] = noisevalue;
			    sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
			}
		    }
		    else {
			if (signalvalue > *cut) {
			    if (caught[i]==0) {              /* first capture of this animal */
				nc++;
				caught[i] = nc;
				for (j=0; j<*ss; j++)
				    for (l=0; l<*kk; l++)
					value[*ss * ((nc-1) * *kk + l) + j] = 0;
			    }
			    nd++;
			    if (nd > maxdet) {
				*resultcode = 2;
				return;  /* error */
			    }
			    value[*ss * ((caught[i]-1) * *kk + k) + s] = 1;
			    worksignal[nd-1] = signalvalue;
			    sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
			}
		    }
                }
            }
        }
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++)
        signal[i]   = worksignal[sortorder[i]];
    if ((*fn == 12) || (*fn == 13)) {
	for (i=0; i<nd; i++)
	    noise[i]   = worknoise[sortorder[i]];
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/


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
    int    *exactn,  /* 0 or a positive integer for the exact number of fixes per animal */
    int    *maxperpoly, /*   */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
)

{
    int    i,j,s,t;
    int    nc = 0;
    int    np = 1;       /* number of points each call of gxy */
    int    nd = 0;
    int    count;
    double par[3];
    double ws;
    int    maxdet;
    double xy[2];
    double *workXY;
    int    *sortorder;
    double *sortkey;
    double Tski = 1.0;   /* not defined !! 2012-12-18 */

    *resultcode = 1;
    GetRNGstate();
    maxdet = *N * *ss * *maxperpoly;
    workXY = (double*) R_alloc(maxdet*2, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));

    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        ws = 10 * sigma[s];
        par[0] = lambda[s];
        par[1] = sigma[s];
        par[2] = z[s];
        if (lambda[s]>0) {
            for (i=0; i<*N; i++) {
		if (*exactn)
		    count = *exactn;
		else
		    count = rcount (binomN[s], lambda[s], Tski); 
                /* require maximum at r=0 */
                if (*fn == 6) error ("annular normal not allowed in trappingtelemetry");
                par[0] = 1;
                for (j=0; j<count; j++) {
                    gxy (&np, fn, par, &ws, xy);            /* simulate location */
                    xy[0] = xy[0] + animals[i];
                    xy[1] = xy[1] + animals[*N + i];
		    if (caught[i]==0) {            /* first capture of this animal */
			nc++;
			caught[i] = nc;
			for (t=0; t<*ss; t++)
			    value[*ss * (nc-1) + t] = 0;
		    }
		    nd++;
		    if (nd >= maxdet) {
			*resultcode = 2;           /* error */
			return;
		    }
		    value[*ss * (caught[i]-1) + s]++;
		    workXY[(nd-1)*2] = xy[0];
		    workXY[(nd-1)*2+1] = xy[1];
		    sortkey[nd-1] = (double) (s * *N + caught[i]);
		}
	    }
	}
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++) {
        detectedXY[i]    = workXY[sortorder[i]*2];
        detectedXY[i+nd] = workXY[sortorder[i]*2+1];
    }
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

void trappingtimes (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    double *dist2,     /* distances squared (optional: -1 if unused) */
    double *Tsk,       /* ss x kk array of 0/1 usage codes or effort */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,        /* truncation radius */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *times,     /* time of detection within occasion */
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
)
{
    double lambda, timevalue;
    int    i,j,k,l,s;
    int    nc = 0;
    int    nd = 0;
    double d2val;
    int    maxdet;
    double *work;
    int    *sortorder;
    double *sortkey;
    double miscparm[3];
    double Tski;

    *resultcode = 1;
    GetRNGstate();
    maxdet = *N * *ss * *kk * 100;
    work = (double*) R_alloc(maxdet, sizeof(double));
    sortorder = (int*) R_alloc(maxdet, sizeof(int));
    sortkey = (double*) R_alloc(maxdet, sizeof(double));
    for (i=0; i<*N; i++) caught[i] = 0;
    for (s=0; s<*ss; s++) {
        for (i=0; i<*N; i++) {
            for (k=0; k<*kk; k++) {
                /* if (used[s * *kk + k]) { */
		Tski = Tsk[s * *kk + k];
		if (fabs(Tski) > 1e-10) {          /* 2012 12 18 */
                    timevalue = 0;
		    d2val = d2L(k, i, dist2, *kk);
                    lambda = pfn(*fn, d2val, g0[k * *ss + s], sigma[s], z[s], miscparm, *w2); 
                    if (lambda>0) {
                        while (timevalue < 1) {
                            timevalue += rexp(1/lambda);
                            if (timevalue < 1) {
                                if (caught[i]==0) {     /* first capture of this animal */
                                    nc++;
                                    caught[i] = nc;
                                    for (j=0; j<*ss; j++)
                                      for (l=0; l<*kk; l++)
                                        value[*ss * ((nc-1) * *kk + l) + j] = 0;
                                }
                                nd++;
                                if (nd >= maxdet) {
                                    *resultcode = 2;
                                    return;  /* error */
                                }
                                value[*ss * ((caught[i]-1) * *kk + k) + s]++;
                                work[nd-1] = timevalue;
                                sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                            }
                        }
                    }
                }
            }
        }
    }
    for (i=0; i<nd; i++)
        sortorder[i] = i;
    if (nd>0) rsort_with_index (sortkey, sortorder, nd);
    for (i=0; i<nd; i++)
        times[i]   = work[sortorder[i]];
    *n = nc;
    *resultcode = 0;
    PutRNGstate();
}
/*==============================================================================*/

