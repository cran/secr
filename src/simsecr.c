/*
    Simulate capture histories from fitted model

    2012-02-13 Tentatively extended for 'signalnoise' detectors 

*/

#include "secr.h"

/*==============================================================================*/

void getpar (int i, int s, int k, int xi, int N, int ss, int nk, int cc0, int cc1, 
             int fn, int bswitch, int gsb0[], double gsb0val[], int gsb1[], 
             double gsb1val[], double *g0, double *sigma, double *z) {
    /* bswitch determines whether to use naive (0) or caught before (1) */
    int wxi;
    int c;
    wxi = i4(i,s,k,xi,N,ss,nk);
    if (bswitch == 0) {
        c = gsb0[wxi]-1;
        *g0 = gsb0val[c];
        *sigma = gsb0val[cc0 + c];
        if ((fn==1) || (fn == 5)  || (fn == 6)  || (fn == 7) || (fn == 8) )
            *z = gsb0val[2* cc0 + c];
    }
    else {
        c = gsb1[wxi]-1;
        *g0 = gsb1val[c];
        *sigma = gsb1val[cc1 + c];
        if ((fn==1) || (fn == 5) || (fn == 6)  || (fn == 7) || (fn == 8) ) 
            *z = gsb1val[2* cc1 + c];
    }
}
/*==============================================================================*/

int rdiscrete (int n, double pmix[])
/* return random discrete observation from distribution in pmix */
{
    double *cumpmix;
    int x;
    double r;
    if (n<1) error ("invalid n in rdiscrete");
    if (n==1) return (0);
    else {
        cumpmix = (double*) R_alloc(n+1, sizeof(double));
        cumpmix[0] = 0;
        for (x=0; x<n; x++) {
            cumpmix[x+1] = cumpmix[x] + pmix[x];
        }
        r = Random();
        for (x=1; x<=n; x++) if (r<cumpmix[x]) break;
        return(x);
    }
}
/*==============================================================================*/

int bswitch (int btype, int N, int i, int k, int caughtbefore[])
{
    if (btype == 0)
        return(0);
    else if (btype == 1) 
        return(caughtbefore[i]);
    else if (btype == 2) 
        return(caughtbefore[k * (N-1) + i]);
    else if (btype == 3) 
        return(caughtbefore[k]);
    else 
        error("unrecognised btype in simsecr");
    return(0);
}
/*==============================================================================*/

void simsecr (
    int    *detect,     /* detector 0 multi, 1 proximity, 2 single, 3 count, 4 area ??? */
    double *gsb0val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    double *gsb1val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [caught before] */
    int    *cc0,        /* number of g0/sigma/b combinations for naive animals */
    int    *cc1,        /* number of g0/sigma/b combinations for caught before */
    int    *gsb0,       /* lookup which g0/sigma/b combination to use for given g, S, K
                           [naive animal] */
    int    *gsb1,       /* lookup which g0/sigma/b combination to use for given n, S, K 
                           [caught before] */
    int    *N,          /* number of animals */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *nmix,       /* number of classes */
    double *animals,    /* x,y points of animal range centres (first x, then y) */
    double *traps,      /* x,y locations of traps (first x, then y) */
    int    *used,       /* ss x kk array of 0/1 codes for usage */
    int    *btype,      /* code for behavioural response 
                           0 none
                           1 individual
                           2 individual, trap-specific
                           3 trap-specific
                        */
    int    *Markov,     /* learned vs transient behavioural response 
                           0 learned
                           1 Markov 
                        */
    int    *binomN,     /* number of trials for 'count' detector modelled with binomial */
    double *miscparm,   /* detection threshold on transformed scale, etc. */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 = uniform */
    int    *maxperpoly, /*  */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* sequence number in session (0 if not caught) */
    double *detectedXY, /* x,y locations of detections  */
    double *signal,     /* vector of signal strengths, one per detection */
    int    *value,      /* return value array of trap locations n x s */
    int    *resultcode
)
{

    double d2val;
    double p;
    int    i,j,k,l,s;
    int    ik;
    int    nc = 0;
    int    nk = 0;             /* number of detectors (polygons or transects when *detect==6,7) */
    int    count = 0;
    int    *caughtbefore;
    int    *x;                 /* mixture class of animal i */
    double *pmix;
    double runif;
    int    wxi = 0;
    int    c = 0;
    int    gpar = 2;
    double g0 = 0;
    double sigma = 0;
    double z = 0;
    double *work = NULL;
    double *noise = NULL;   /* detectfn 12,13 only */
    int    *sortorder = NULL;
    double *sortkey = NULL;

    /*
        *detect may take values -
       -1  single-catch traps
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        5  signal detectors
        6  polygon detectors
        7  transect detectors
    */
    /*========================================================*/
    /* 'single-catch only' declarations */
    int    tr_an_indx = 0;
    int    nanimals;
    int    ntraps;
    int    *occupied = NULL;
    int    *intrap = NULL;
    struct trap_animal *tran = NULL;
    double event_time;
    int    anum = 0;
    int    tnum = 0;
    int    nextcombo;
    int    finished;
    int    OK;

    /*========================================================*/
    /* 'multi-catch only' declarations */
    double *h = NULL;          /* multi-catch only */
    double *hsum = NULL;       /* multi-catch only */
    double *cump = NULL;       /* multi-catch only */

    /*========================================================*/
    /* 'polygon & transect only' declarations */
    int    nd = 0;
    int    cumk[maxnpoly+1];
    int    sumk;               /* total number of vertices */
    int    g=0;
    int    *gotcha;
    double xy[2];
    int    n1,n2,t;
    double par[3];
    int    np = 1;             /* n points each call of gxy */
    double w, ws;
    int    maxdet=1;
    double *cumd = NULL;
    struct rpoint *line = NULL;
    struct rpoint xyp;
    struct rpoint animal;
    double lx;
    double maxg = 0;
    double lambdak;  /* temp value for Poisson rate */
    double grx;      /* temp value for integral gr */
    double stdint;
    int    J;
    int    maybecaught;
    double dx,dy,d;
    double pks;
    double sumhaz;
    /*========================================================*/
    /* 'signal-strength only' declarations */
    double beta0;
    double beta1;
    double muS;
    double sdS;
    double muN = 0;
    double sdN = 1;
    double signalvalue;
    double noisevalue;
    double cut;

    /*========================================================*/
    /* MAIN LINE */

    gotcha = &g;
    *resultcode = 1;
    caughtbefore = (int *) R_alloc(*N * *kk, sizeof(int));
    x = (int *) R_alloc(*N, sizeof(int));
    for (i=0; i<*N; i++) x[i] = 0;
    pmix = (double *) R_alloc(*nmix, sizeof(double));

    if ((*detect < -1) || (*detect > 7)) return;

    if (*detect == -1) {                                   /* single-catch only */
        occupied = (int*) R_alloc(*kk, sizeof(int));
        intrap = (int*) R_alloc(*N, sizeof(int));
        tran = (struct trap_animal *) R_alloc(*N * *kk,
            sizeof(struct trap_animal));
	    /*  2*sizeof(int) + sizeof(double)); */
    }
    if (*detect == 0) {                                    /* multi-catch only */
        h = (double *) R_alloc(*N * *kk, sizeof(double));
        hsum = (double *) R_alloc(*N, sizeof(double));
        cump = (double *) R_alloc(*kk+1, sizeof(double));
        cump[0] = 0;
    }

    if (*detect == 5) {                                    /* signal only */
        maxdet = *N * *ss * *kk;
        if (!((*fn == 10) || (*fn == 11)))
            error ("simsecr not implemented for this combination of detector & detectfn");

    }

    if ((*detect == 3) || (*detect == 4) || (*detect == 6) || (*detect == 7)) {
        /* polygon or transect */
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {                       /* maxnpoly much larger than npoly */
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
        sumk = cumk[nk];
        if ((*detect == 6) || (*detect == 7))
            maxdet = *N * *ss * nk * *maxperpoly;
        else 
            maxdet = *N * *ss;
    }
    else nk = *kk;

    if ((*detect == 4) || (*detect == 7)) {                            /* transect only */
        line = (struct rpoint *) R_alloc(sumk, sizeof(struct rpoint));
        cumd = (double *) R_alloc(sumk, sizeof(double));
        /* coordinates of vertices */
        for (i=0; i<sumk; i++) {
            line[i].x = traps[i];
            line[i].y = traps[i+sumk];
        }
        /* cumulative distance along line; all transects end on end */
        for (k=0; k<nk; k++) {
            cumd[cumk[k]] = 0;
            for (i=cumk[k]; i<(cumk[k+1]-1); i++) {
                cumd[i+1] = cumd[i] + distance(line[i], line[i+1]);
            }
        }
    }

    if ((*detect==3) || (*detect==4) || (*detect==5) || (*detect==6) || (*detect==7)) {
        work = (double*) R_alloc(maxdet*2, sizeof(double));   /* twice size needed for signal */
        sortorder = (int*) R_alloc(maxdet, sizeof(int));
        sortkey = (double*) R_alloc(maxdet, sizeof(double));
    }
    if ((*fn==12) || (*fn==13)) {
        noise = (double*) R_alloc(maxdet*2, sizeof(double));   /* twice size needed for signal */
    }

    GetRNGstate();

    /* may be better to pass pmix */
    gpar = 2;
    if ((*fn == 1) || (*fn == 3) || (*fn == 5)|| (*fn == 6) || (*fn == 7) || (*fn == 8) ||
        (*fn == 10) || (*fn == 11)) gpar ++;
    if (*nmix>1) gpar++;

    if (*nmix>1) {
        if (*nmix>2)
            error("simsecr nmix>2 not implemented");
        for (i=0; i<*nmix; i++) {
            wxi = i4(0,0,0,i,*N,*ss,nk);
            c = gsb0[wxi] - 1;
            pmix[i] = gsb0val[*cc0 * (gpar-1) + c];    /* assuming 4-column gsb */
        }
        for (i=0; i<*N; i++) {
            x[i] = rdiscrete(*nmix, pmix) - 1;
        }
    }

    /* zero caught status */
    for (i=0; i<*N; i++) 
        caught[i] = 0;    
    for (i=0; i<*N; i++) 
        for (k=0; k < nk; k++)
            caughtbefore[k * (*N-1) + i] = 0;

    /* MAIN LOOP */
    for (s=0; s<*ss; s++) {

        /* ------------------ */
        /* single-catch traps */
        if (*detect == -1) {
            /* initialise day */
            tr_an_indx = 0;
            nanimals = *N;
            ntraps   = nk;
            for (i=0; i<*N; i++) intrap[i] = 0;
            for (k=0; k<nk; k++) occupied[k] = 0;
            nextcombo = 0;

            /* make tran */
            for (i=0; i<*N; i++) {  /* animals */
                for (k=0; k<nk; k++) { /* traps */
                    if (used[s * nk + k]) {                        
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
			    bswitch (*btype, *N, i, k, caughtbefore), 
                            gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
                        d2val = d2(i,k, animals, traps, *N, nk);
                        p = pfn(*fn, d2val, g0, sigma, z, miscparm, 1e20);   /* effectively infinite w2 */
                        event_time = randomtime(p);
                        if (event_time <= 1) {
                            tran[tr_an_indx].time   = event_time;
                            tran[tr_an_indx].animal = i;    /* 0..*N-1 */
                            tran[tr_an_indx].trap   = k;    /* 0..nk-1 */
                            tr_an_indx++;
                        }
                    }
                }
            }
            /* end of make tran */

            if (tr_an_indx > 1) probsort (tr_an_indx, tran);

            while ((nextcombo < tr_an_indx) && (nanimals>0) && (ntraps>0)) {
                    finished = 0;
                OK       = 0;
                    while ((1-finished)*(1-OK) > 0) {      /* until finished or OK */
                    if (nextcombo >= (tr_an_indx))
                            finished = 1;                  /* no more to process */
                    else {
                            anum = tran[nextcombo].animal;
                        tnum = tran[nextcombo].trap;
                            OK = (1-occupied[tnum]) * (1-intrap[anum]); /* not occupied and not intrap */
                        nextcombo++;
                        }
                }
                    if (finished==0) {
                       /* Record this capture */
                          occupied[tnum] = 1;
                      intrap[anum]   = tnum+1;         /* trap = k+1 */
                          nanimals--;
                      ntraps--;
                    }
            }
                for (i=0; i<*N; i++) {
                if (intrap[i]>0) {
                    if (caught[i]==0) {                    /* first capture of this animal */
                       nc++;
                       caught[i] = nc;                     /* nc-th animal to be captured */
                       for (j=0; j<*ss; j++)
                           value[*ss * (nc-1) + j] = 0;
                    }
                    value[*ss * (caught[i]-1) + s] = intrap[i];  /* trap = k+1 */
                }
            }
        }

        /* -------------------------------------------------------------------------- */
        /* multi-catch trap; only one site per occasion (drop last dimension of capt) */
        else if (*detect == 0) {
            for (i=0; i<*N; i++) {
                hsum[i] = 0;
                for (k=0; k<nk; k++) {
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
			    bswitch (*btype, *N, i, k, caughtbefore), 
                            gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
                    d2val = d2(i,k, animals, traps, *N, nk);
                    p = pfn(*fn, d2val, g0, sigma, z, miscparm, 1e20);
                    p = p * used[s * nk + k];           /* zero if not used */
                    h[k * *N + i] = -log(1 - p);
                    hsum[i] += h[k * *N + i];
                }

                for (k=0; k<nk; k++) {
                    cump[k+1] = cump[k] + h[k * *N + i]/hsum[i];
                }
                if (Random() < (1-exp(-hsum[i]))) {
                    if (caught[i]==0)  {        /* first capture of this animal */
                        nc++;
                        caught[i] = nc;
                        for (j=0; j<*ss; j++)
                            value[*ss * (nc-1) + j] = 0;
                    }
                    /* find trap with probability proportional to p
                       searches cumulative distribution of p  */
                    runif = Random();
                    k = 0;
                    while ((runif > cump[k]) && (k<nk)) k++;
                    value[*ss * (caught[i]-1) + s] = k;  /* trap = k+1 */
                }
            }
        }

        /* -------------------------------------------------------------------------------- */
        /* the 'proximity' group of detectors 1:2 - proximity, count */
        else if ((*detect >= 1) && (*detect <= 2)) {
            for (i=0; i<*N; i++) {
                for (k=0; k<nk; k++) {
                    if (used[s * nk + k]) {
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
			    bswitch (*btype, *N, i, k, caughtbefore), 
                            gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
                        d2val = d2(i,k, animals, traps, *N, nk);
                        p = pfn(*fn, d2val, g0, sigma, z, miscparm, 1e20);

                        if (p < -0.1) { PutRNGstate(); return; }   /* error */

                        if (p>0) {
                            if (*detect == 1) {
                                count = Random() < p;              /* binary proximity */
                            }
                            else if (*detect == 2) {               /* count proximity */
                                count = rcount(*binomN, p);
                            }
                            if (count>0) {
                                if (caught[i]==0) {                /* first capture of this animal */
                                    nc++;
                                    caught[i] = nc;
                                    for (j=0; j<*ss; j++)
                                      for (l=0; l<nk; l++)
                                        value[*ss * ((nc-1) * nk + l) + j] = 0;
                                }
                                value[*ss * ((caught[i]-1) * nk + k) + s] = count;
                            }
                        }
                    }
                }
            }
        }

        /* -------------------------------------------------------------------------------- */
        /* exclusive polygon detectors  */
        else if (*detect == 3) {
            /* find maximum distance between animal and detector vertex */
            w = 0;
            J = cumk[nk];
            for (i = 0; i< *N; i++) {
                for (j = 0; j < J; j++) {
                    dx = animals[i] - traps[j];
                    dy = animals[*N + i] - traps[J + j];
        	    d = sqrt(dx*dx + dy*dy);
                    if (d > w) w = d;
                }
            } 

            for (i=0; i<*N; i++) {
                /* this implementation assumes NO VARIATION AMONG DETECTORS */
                getpar (i, s, 0, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
		    bswitch (*btype, *N, i, 0, caughtbefore), 
                    gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
                maybecaught = Random() < g0;

                if (w > (10 * sigma)) 
                    ws = 10 * sigma;
                else 
                    ws = w;

                par[0] = 1;
                par[1] = sigma;
                par[2] = z;
                
                if (maybecaught) {
/* 2011-05-10          gxy (&np, fn, par, &w, xy);   */
                    gxy (&np, fn, par, &ws, xy);                 /* simulate location */
                    xy[0] = xy[0] + animals[i];
                    xy[1] = xy[1] + animals[*N + i];
                    for (k=0; k<nk; k++) {                      /* each polygon */
                        if (used[s * nk + k]) {
                            n1 = cumk[k];
                            n2 = cumk[k+1]-1;
                            inside(xy, &n1, &n2, &sumk, traps, gotcha);  /* assume closed */
                            if (*gotcha > 0) {
                                if (caught[i]==0) {             /* first capture of this animal */
                                    nc++;
                                    caught[i] = nc;
                                    for (t=0; t<*ss; t++)
                                            value[*ss * (nc-1) + t] = 0;
                                }
                                nd++;
                                value[*ss * (caught[i]-1) + s] = k+1;
                                work[(nd-1)*2] = xy[0];
                                work[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (s * *N + caught[i]);
                                break;   /* no need to look at more poly */
                            }
                        }
                    }
                }
            }
        }

        /* -------------------------------------------------------------------------------- */
        /* exclusive transect detectors  */
        else if (*detect == 4) {
            for (i=0; i<*N; i++) {                            /* each animal */
                animal.x = animals[i];
                animal.y = animals[i + *N];
                sumhaz = 0;
                /* ------------------------------------ */
                /* sum hazard */
                for (k=0; k<nk; k++) {            
                    if (used[s * nk + k]) {
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
                            bswitch (*btype, *N, i, k, caughtbefore), 
                            gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
	                par[0] = g0;
                        par[1] = sigma;
                        par[2] = z;
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        stdint = gintegral1(*fn, par);
                        sumhaz += -log(1 - par[0] * integral1D (*fn, i, 0, par, 1, traps, 
                            animals, n1, n2, sumk, *N) / stdint);
                    }
                }
                /* ------------------------------------ */

                for (k=0; k<nk; k++) {                        /* each transect */
                    if (used[s * nk + k]) {
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
   			    bswitch (*btype, *N, i, k, caughtbefore), 
                            gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
	                par[0] = g0;
                        par[1] = sigma;
                        par[2] = z;
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        stdint = gintegral1(*fn, par);
                        lambdak = par[0] * integral1D (*fn, i, 0, par, 1, traps, animals,
                            n1, n2, sumk, *N) / stdint;
    	                pks = (1 - exp(-sumhaz)) * (-log(1-lambdak)) / sumhaz;
                        count = Random() < pks;
                        maxg = 0;
                        if (count>0) {                       /* find maximum - approximate */
                            for (l=0; l<=100; l++) {
                                lx = (cumd[n2] - cumd[n1]) * l/100;
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = fmax2(maxg, grx);
                            }
                            for (l=n1; l<=n2; l++) {
                                xyp = line[l];
                                grx = gr (fn, par, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = fmax2(maxg, grx);
                            }
                            maxg= 1.2 * maxg;                 /* safety margin */
                            if (maxg<=0)
                                Rprintf("maxg error in simsecr\n"); /* not found */
                            *gotcha = 0;
                            l = 0;
                            while (*gotcha == 0) {
                                lx = Random() * (cumd[n2] - cumd[n1]);     /* simulate location */
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xyp, animal);
                                if (Random() < (grx/maxg))    /* rejection sampling */
                                    *gotcha = 1;
                                l++;
                                if (l % 10000 == 0)
                                    R_CheckUserInterrupt();
                                if (l>1e8) *gotcha = 1;       /* give up and accept anything!!!! */
                            }
                            if (caught[i]==0) {               /* first capture of this animal */
                                nc++;
                                caught[i] = nc;
                                for (t=0; t<*ss; t++)
                                        value[*ss * (nc-1) + t] = 0;
                            }
                            nd++;
                            if (nd >= maxdet) {
                                *resultcode = 2;  /* error */
                                return;
                            }
                            value[*ss * (caught[i]-1) + s] = k+1;
                            work[(nd-1)*2] = xyp.x;
                            work[(nd-1)*2+1] = xyp.y;
                            sortkey[nd-1] = (double) (s * *N + caught[i]);
                        }
                        if (count>0) break;   /* no need to look further */
                    }
                }                                             /* end loop over transects */
            }                                                 /* end loop over animals */
        }

        /* -------------------------------------------------------------------------------- */
        /* polygon detectors  */
        else if (*detect == 6) {
            for (i=0; i<*N; i++) {
                /* this implementation assumes NO VARIATION AMONG DETECTORS */

                getpar (i, s, 0, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
		    bswitch (*btype, *N, i, 0, caughtbefore), 
                    gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
                count = rcount(*binomN, g0);
                w = 10 * sigma;
                par[0] = 1;
                par[1] = sigma;
                par[2] = z;
                for (j=0; j<count; j++) {
                    gxy (&np, fn, par, &w, xy);                 /* simulate location */
                    xy[0] = xy[0] + animals[i];
                    xy[1] = xy[1] + animals[*N + i];
                    for (k=0; k<nk; k++) {                      /* each polygon */
                        if (used[s * nk + k]) {
                            n1 = cumk[k];
                            n2 = cumk[k+1]-1;
                            inside(xy, &n1, &n2, &sumk, traps, gotcha);  /* assume closed */
                            if (*gotcha > 0) {
                                if (caught[i]==0) {             /* first capture of this animal */
                                    nc++;
                                    caught[i] = nc;
                                    for (t=0; t<*ss; t++)
                                        for (l=0; l<nk; l++)
                                            value[*ss * ((nc-1) * nk + l) + t] = 0;
                                }
                                nd++;
                                if (nd > maxdet) {
                                    *resultcode = 2;
                                    return;  /* error */
                                }
                                value[*ss * ((caught[i]-1) * nk + k) + s]++;
                                work[(nd-1)*2] = xy[0];
                                work[(nd-1)*2+1] = xy[1];
                                sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                            }
                        }
                    }
                }
            }
        }

        /* -------------------------------------------------------------------------------- */
        /* transect detectors  */
        else if (*detect == 7) {
            for (i=0; i<*N; i++) {                            /* each animal */
                animal.x = animals[i];
                animal.y = animals[i + *N];
                for (k=0; k<nk; k++) {                        /* each transect */
                    if (used[s * nk + k]) {
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
 			    bswitch (*btype, *N, i, k, caughtbefore), 
                            gsb0, gsb0val, gsb1, gsb1val, &g0, &sigma, &z);
	                par[0] = g0;
                        par[1] = sigma;
                        par[2] = z;
                        n1 = cumk[k];
                        n2 = cumk[k+1]-1;
                        stdint = gintegral1(*fn, par);
                        lambdak = par[0] * integral1D (*fn, i, 0, par, 1, traps, animals,
                            n1, n2, sumk, *N) / stdint;
                        count = rcount(*binomN, lambdak);    /* number of detections on transect */
                        maxg = 0;
                        if (count>0) {                       /* find maximum - approximate */
                            for (l=0; l<=100; l++) {
                                lx = (cumd[n2]-cumd[n1]) * l/100;
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = fmax2(maxg, grx);
                            }
                            for (l=n1; l<=n2; l++) {
                                xyp = line[l];
                                grx = gr (fn, par, xyp, animal);
                                if (R_FINITE(grx))
                                    maxg = fmax2(maxg, grx);
                            }
                            maxg= 1.2 * maxg;                 /* safety margin */
                            if (maxg<=0)
                                Rprintf("maxg error in simsecr\n"); /* not found */

                        }
                        for (j=0; j<count; j++) {
                            *gotcha = 0;
                            l = 0;
                            while (*gotcha == 0) {
                                lx = Random() * (cumd[n2]-cumd[n1]);     /* simulate location */
                                xyp = getxy (lx, cumd, line, sumk, n1);
                                grx = gr (fn, par, xyp, animal);
                                if (Random() < (grx/maxg))   /* rejection sampling */
                                    *gotcha = 1;
                                l++;
                                if (l % 10000 == 0)
                                    R_CheckUserInterrupt();
                                if (l>1e8) *gotcha = 1;        /* give up and accept anything!!!! */
                            }
                            if (caught[i]==0) {               /* first capture of this animal */
                                nc++;
                                caught[i] = nc;
                                for (t=0; t<*ss; t++)
                                    for (l=0; l<nk; l++)
                                        value[*ss * ((nc-1) * nk + l) + t] = 0;
                            }
                            nd++;
                            if (nd >= maxdet) {
                                *resultcode = 2;  /* error */
                                return;
                            }
                            value[*ss * ((caught[i]-1) * nk + k) + s]++;
                            work[(nd-1)*2] = xyp.x;
                            work[(nd-1)*2+1] = xyp.y;
                            sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
                        }
                    }
                }                                             /* end loop over transects */
            }                                                 /* end loop over animals */
        }

        /* ------------------------ */
        /* signal strength detector */
        else if (*detect == 5) {
	    cut = miscparm[0];
	    if ((*fn == 12) || (*fn == 13)) {
		muN = miscparm[1];
		sdN = miscparm[2];
	    }
            for (i=0; i<*N; i++) {
                for (k=0; k<nk; k++) {
                    if (used[s * nk + k]) {
                        /* sounds not recaptured */
                        getpar (i, s, k, x[i], *N, *ss, nk, *cc0, *cc1, *fn, 
				0, gsb0, gsb0val, gsb0, gsb0val, &beta0, &beta1, &sdS);
                        if ((*fn == 10) || (*fn == 12))
			    muS  = mufn (i, k, beta0, beta1, animals, traps, *N, nk, 0);
                        else
                            muS  = mufn (i, k, beta0, beta1, animals, traps, *N, nk, 1);
			signalvalue = norm_rand() * sdS + muS;
                        if ((*fn == 10) || (*fn == 11)) {
			    if (signalvalue > cut) {
				if (caught[i]==0) {        /* first capture of this animal */
				    nc++;
				    caught[i] = nc;
				    for (j=0; j<*ss; j++)
					for (l=0; l<nk; l++)
					    value[*ss * ((nc-1) * *kk + l) + j] = 0;
				}
				nd++;
				value[*ss * ((caught[i]-1) * *kk + k) + s] = 1;
				work[nd-1] = signalvalue;
				sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
			    }
			}
			else {
			    noisevalue = norm_rand() * sdN + muN;
			    if ((signalvalue - noisevalue) > cut) {
				if (caught[i]==0) {        /* first capture of this animal */
				    nc++;
				    caught[i] = nc;
				    for (j=0; j<*ss; j++)
					for (l=0; l<nk; l++)
					    value[*ss * ((nc-1) * *kk + l) + j] = 0;
				}
				nd++;
				value[*ss * ((caught[i]-1) * *kk + k) + s] = 1;
				work[nd-1] = signalvalue;
				noise[nd-1] = noisevalue;
				sortkey[nd-1] = (double) (k * *N * *ss + s * *N + caught[i]);
			    }
			}
		    }
                }
            }
        }

	if ((*btype > 0) && (s < (*ss-1))) {
            /* update record of 'previous-capture' status */
            if (*btype == 1) {
                for (i=0; i<*N; i++) {
                    if (*Markov) 
                        caughtbefore[i] = 0;
                    for (k=0; k<nk; k++)
                        caughtbefore[i] = imax2 (value[i3(s, k, i, *ss, nk)], caughtbefore[i]);
                }
            }
            else if (*btype == 2) {
                for (i=0; i<*N; i++) {
                    for (k=0; k<nk; k++) {
                        ik = k * (*N-1) + i;
                        if (*Markov) 
                            caughtbefore[ik] = value[i3(s, k, i, *ss, nk)];
                        else 
                            caughtbefore[ik] = imax2 (value[i3(s, k, i, *ss, nk)], 
			        caughtbefore[ik]);
		    }
		}
            }
	    else {
                for (k=0;k<nk;k++) {
                    if (*Markov) 
                        caughtbefore[k] = 0;
                    for (i=0; i<*N; i++) 
                        caughtbefore[k] = imax2 (value[i3(s, k, i, *ss, nk)], caughtbefore[k]);
                }
	    }
	}

    }   /* loop over s */

    if ((*detect==3) || (*detect==4) || (*detect==5) || (*detect==6) || (*detect==7)) {
        for (i=0; i<nd; i++) sortorder[i] = i;
        if (nd>0) rsort_with_index (sortkey, sortorder, nd);
        if (*detect==5) {
            for (i=0; i<nd; i++) signal[i] = work[sortorder[i]];
	    if ((*fn == 12) || (*fn == 13)) {
		for (i=0; i<nd; i++) signal[i+nd] = noise[sortorder[i]];
	    }
	}
        else {
            for (i=0; i<nd; i++) {
                detectedXY[i]    = work[sortorder[i]*2];
                detectedXY[i+nd] = work[sortorder[i]*2+1];
            }
        }
    }
    *n = nc;
    PutRNGstate();
    *resultcode = 0;
}

/*==============================================================================*/

