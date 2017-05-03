/*
   External procedures for secr package

   can compile with gcc 4.6.3 :
   gcc -Ic:/R/R-3.3.0/include -c secr.c -Wall -pedantic -std=gnu99

   [confirmed 2014-08-27 1021]

*/
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
/* 2012-11-13 moved some code to utils.c and polygon.c */
/* 2012-11-13 code shared by integralprw1, secrloglik and pwuniform tidied */
/* 2012-11-24 fixed polygon bug from hangover switch0 line in precompute */
/* 2012-12-17 Tsk effort adjustment */
/* 2012-12-23 radical hack - remove switch0 etc etc */
/* 2012-12-25 cleaned */
/* 2013-04-14 known classes (hcov) for CL = T */
/* 2013-04-16 change group indexing to zero-based */
/* 2013-04-17 known classes (hcov) for CL = F (incompatible with groups) */
/* 2013-06-05 known class coded 1 for unknown, >1 for known */
/* 2013-07-02 slight adjustments for param=2 (treated as param=0) Dropped 2015-11-19 */
/* 2013-11-16 purged 'nested' cuerate code */
/* 2013-11-23 pimask used for individual distribution of location if pimask[0] >-tol */
/* 2013-12-01 like = 4 for CL concurrent telemetry; anyvarying funtion */
/* 2014-04-03 normalization with getdenom not implemented because of memory issue  */
/*            see test.R, valgrind etc. */
/* 2014-08-27 major rejig to use lookup distances (d2L) computed with makedist2 */
/*            and saved in dist2. This is preliminary to user-provided distance matrix */
/*            These functions do NOT use d2L and dist2 */
/*            -- prwi fns for polygons and transects, and pdotpoly */
/* 2014-09-10 fxIHP drastically redone - now much cleaner */
/* 2014-10-04 geth drastic speed improvement for constant detector covariates WITHDRAWN */
/* 2015-10-01 pndotgrp cuerate function removed */
/* 2015-10-04 updated arguments of prwfn's to accommodate MR; dropped s1,s2*/
/* 2015-10-10 bug fixed: conflict h2 and 3-parameter detection functions */
/* 2015-10-11 fxIHP mark-resight working */
/* 2015-10-22 like = 5,6 for sighting-only */
/* 2015-11-19 dropped param = 1 GR */
/* 2015-11-20 fixed bug in fxIHP that crashed R when detector type='multi' and usage varied */
/* 2015-12-21 temporary fix prwicount */
/* 2016-01-08 correct overflow and make memrequest in secrloglik more specific */
/* 2016-10-06 secr 3.0 -- detect as an integer vector */
/* 2016-10-07 secr 3.0 -- remove Gardner et al 2009 param fn = 20 */
/* 2016-12-30 remove superbinomial option distrib = 2 */
/* 2017-01-05 remove like = 3,4 (CT) and usepimask */
/* 2017-03-23 stdint renamed H */
/* 2017-03-23 gintegral renamed hintegral */
/* 2017-04-03 fixed bug for polygons in getdetspec */
/* 2017-04-03 only precompute gk0 if gsb0val differs from gsbval (behavioural effects) */
/* 2017-04-04 hdotpoly replaces pdotpoly */
/*
        detect[s] may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  exclusive polygon detector
        4  exclusive transect detector
        5  signal detector
        6  polygon detector
        7  transect detector
	10 unmarked
	11 presence/absence
     	12 signalnoise
	13 telemetry
*/

/* *like codes for the type of likelihood. It may take values -

0 full
1 conditional
2 full, split integral for sightings of marked animals (from 2016-11-23)
5 all-sighting histories, n0 known (includes all-zero histories)
6 all-sighting histories, n0 unknown

*/

/* 
2016-12-30 Development of likelihood component for telemetry fixes 

Restricted to detectfn 14,16 (HHN, HEX)
Supercedes previous clunky CT and DT implementation in like=3, like=4

Function getrhr -- precomputes cc x mm x nt array of pdf components 'rhr',
                   and the normalising factor 'normal'
                -- called after 'precompute'
Function fnu    -- probability density for the fixes of one animal 
                -- called after prwi

2017-03-24 Cleanup of polygon hazard formulation

Complete reliance on integrating hazard (zfn)
Removed unnecessary factor of 10000 from prwipolygon, prwipolygonX
Removed unnecessary factor of 100 from prwitransect and prwitransectX

 */

#include "secr.h"
#include <time.h>

/*==============================================================================*/

void R_CheckUserInterrupt(void);

/*==============================================================================*/

/* Functions to characterize detector type */
/* polygon, transect and signal detector types must be constant 
   across occasions */

int anyexclusive (int detect[], int ss) {
    int s;
    int exclusive = 0;
    for (s=0; s<ss; s++) 
	if ((detect[s]==0) || (detect[s]==3) || (detect[s]==4))
	    exclusive = 1;
    return exclusive;
}
int anypolygon (int detect[], int ss) {
    int s;
    int polygon = 0;
    for (s=0; s<ss; s++) 
	if ((detect[s]==3) || (detect[s]==6) )
	    polygon = 1;
    return polygon;
}
int anytransect (int detect[], int ss) {
    int s;
    int transect = 0;
    for (s=0; s<ss; s++) 
	if ((detect[s]==4) || (detect[s]==7))
	    transect = 1;
    return transect;
}
int anysignal (int detect[], int ss) {
    int s;
    int signal = 0;
    for (s=0; s<ss; s++) 
	if ((detect[s]==5) || (detect[s]==12))
	    signal = 1;
    return signal;
}
int anytelemetry (int detect[], int ss) {
    int s;
    int telemetry = 0;
    for (s=0; s<ss; s++) 
	if ((detect[s]==13))
	    telemetry = 1;
    return telemetry;
}
int alltelemetry (int detect[], int ss) {
    int s;
    int telemetry = 1;
    for (s=0; s<ss; s++) 
	if ((detect[s]!=13))
	    telemetry = 0;
    return telemetry;
}
int allpoint (int detect[], int ss, int allowsignal, int allowtelem) {
    int s;
    int point;
    int OK = 1;
    for (s=0; s<ss; s++) {
	point = (detect[s]==0) || (detect[s]==1) || (detect[s]==2)
	    || (detect[s]==10) || (detect[s]==11)
	    || (allowsignal && ((detect[s]==5) || (detect[s]==12)))
	    || (allowtelem && ((detect[s]==13) || (detect[s]==8)));
	OK = OK && point;
    }
    return OK;
}
int allmulti (int detect[], int ss) {
    int s;
    int notmulti = 0;
    for (s=0; s<ss; s++) 
	if (detect[s]!=0)
	    notmulti = 1;
    return 1-notmulti;
}

/* Do parameter values for naive animals differ at all from thos for other animals ?*/
int anyb (double gsbval[], double gsb0val[], int cc, int npar) {
    int identical = 1;
    int i;
    for (i=0; i<(cc*npar); i++) {
	if (gsbval[i] != gsb0val[i]) identical = 0;
    }
    return 1 - identical;
}
/*==============================================================================*/

double pndot (int m, int n, int markocc[], int x, int ncol, int PIA0[],
	      double gk0[], double hk0[], int detect[], int binomN[], double Tsk[], int ss, int nk, 
	      int cc0, int nmix, double gsb0val[], int allsighting)

/*
    probability animal at point m on mask is caught
    n may indicate group (full likelihood; ncol= number of groups) or
    individual (conditional likelihood; ncol= number of individuals)
    aligned with secrloglik 2009 06 25

    2009 10 24 adjusted to allow summation over qq < ss
    2009 11 12 'nk' should be number of parts for polygon detectors
    2011 01 04 'param = 1' for GR parameterisation of multi Dropped 2015-11-19
    2012 12 17 Tsk
    2012 12 25 for x=0, dbinom(x, N, 1-(1-p)^Tsk) = dbinom(x, N * Tsk, p)
    2012 12 25 DOES NOT ALLOW binomN < 0
    2015-03-18 used 'p0' as variable for Pr(not caught) instead of re-using g1
    2015-11-19 dropped 'param=1'
    2016-12-31 notional telemetry trap (nk) is always usage=0 for detect[s] != 13
    2017-03-20 hk; force polygon = count
*/
{
    int k,s,c;
    int gi;
    int wxi;
    double pp;
    double p0;
    double Tski;
    pp = 1;
    for (s=0; s<ss; s++) {
        if (((markocc[s]>0) || allsighting) && (detect[s]!=13)) {   /* capture occasions */
	    if (binomN[s] < 0) error ("negative binomN not allowed in C fn pndot");
	    for (k=0; k< nk; k++) {
		wxi = i4(n,s,k,x,ncol,ss,nk);
		c = PIA0[wxi] - 1;
		if (c >= 0) {    /* drops unset traps */
		    gi = i3(c,k,m,cc0,nk);
		    Tski = Tsk[s * nk + k];
		    /* expect binomN = 1 if not count detector */
		    if (fabs(Tski-1) > 1e-10) {                  /* effort <> 1.0 */
			if ((detect[s] < 8) & (detect[s] != 5))  {
			    if (binomN[s] == 0)
				p0 = exp(-Tski *  hk0[gi]);      /* Poisson count or polygon detr */
			    else
				p0 = pow(1-gk0[gi], Tski * binomN[s]);  /* Binomial and Bernoulli */
			}
			else error("no effort adjustment for detector type");
		    }
		    else {
			if (binomN[s] == 0) 
				p0 = exp(- hk0[gi]);       /* Poisson count or polygon detector */
			else 
			    if (binomN[s] == 1)
			    p0 = 1 - gk0[gi];                   /* Bernoulli */
			else 
			    p0 = pow(1-gk0[gi], binomN[s]);     /* Binomial */			
			
		    }
		    pp *= p0;
		    
		}
	    }
	}
    }
    /* Rprintf("pp in pndot %12.6f\n", pp); */
    return (1 - pp);    
}
/*===============================================================*/

double prwipoint
    (int m, int n, int x, int w[], double xy[], double signal[],
     int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
     int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
     double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])
    
    /*

    Probability of capture history n (0 <= n < *nc)
    given that animal's range centre is at m
    
    Combined code for all point detectors (previously prwimulti, prwiprox, prwicount)

    */
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of trap      0 <= k < kk  */
    int c, gi, wi, w0, wxi;
    int count;
    int detect;
    int dead = 0;
    double htemp;
    double result = 1.0;
    double g1;
    double Tski;
    double pks;
    double pI = 1.0;
    for (s=0; s<ss; s++) {
	pI = pID[s + ss*x];
	detect = (int) (detspec[s]+0.1);
	// if (m==0) Rprintf("prwipoint s %3d detect %3d\n", s, detect);
        /* sneaky way to control which occasions to use */
	if ((detect>=0) && (detect!=13)) {   

	    if (detect == 0) {                                     /* trap */
                /* assume trapped animals are never 'sightings'; ignore pI */
		htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
		if (htemp < fuzz) { result = 0; break; }
		w0 = i3(n, s, 0, nc, ss);
		if (w[w0] != 0) {                                /* Captured */
		    // if(m==0) Rprintf("n %5d s %5d  w[w0] %5d  \n", n, s, w[w0]);
		    if (w[w0] < 0) dead = 1;  
		    k = abs(w[w0])-1;
		    Tski = Tsk[s * kk + k];
		    wxi = i4(n, s, k, x, nc, ss, kk);
		    c = PIA[wxi] - 1;
		    gi = i3(c, k, m, cc, kk);
		    pks = Tski * hk[gi];
		    pks *=  (1-expmin(-htemp)) / htemp;
		}
		else  {
		    pks = expmin(-htemp);      /* Not captured */
		}
		result *= pks;	
	    }
	    else {
		for (k=0; k<kk; k++) {
		    wi = i3(n, s, k, nc, ss);
		    wxi = i4(n, s, k, x, nc, ss, kk);
		    count = w[wi];
		    if (count<0) {dead=1; count=-count;}
		    c = PIA[wxi] - 1;
		    if (c >= 0) {                                     /* skip if this trap not set */
			gi  = i3(c, k, m, cc, kk);
			g1 = gk[gi];
			Tski = Tsk[s * kk + k];

			if (fabs(Tski-1) > 1e-10) {                   /* not unity */ 
			    if (Tski < 1e-10)                         /* 2014-11-10 */
				g1 = 0;
			    else
				g1 = 1 - pow(1 - g1, Tski);
			}
			if (detect == 1 ) {                           /* binary proximity    */
			    if (count)                                /* Bernoulli count 0/1 */
				result *= g1 * pI;
			    else 
				result *= (1 - g1 * pI);
			}
			else {                                        /* count proximity */

			    if (binomN[s] == 0) {                     /* Poisson */
				if (count == 0) 
				    result *= expmin(-Tski * hk[gi] * pI);
				else
				    result *= dpois(count, Tski * hk[gi] * pI, 0); 
			    }
			    else if (binomN[s] == 1) {                /* Binomial, size from Tsk */ 
				result *= countp (count, round(Tski), g1 * pI);
			    }
			    else if (binomN[s] > 1) {                 /* Binomial, specified size */
				if (fabs(Tski-1) > 1e-10)             /* not unity */
				    g1 = 1 - pow(1 - g1, Tski);
				result *= countp (count, binomN[s], g1 * pI);
			    }
			}
			if (result < minp) {result = minp; break;}    /* truncate */
		    }
		}
		if (result <= minp) {result = minp; break;}           /* truncate */
	    }
	    if (dead==1) break;
	}
    }
    return (result);
}
/*=============================================================*/

void gethr(int detect[], int fn, int ss, int mm, int cc, int nt,
	    double mask[], double xy[], double gsbval[], double hr[],
            double telemscale) {

    /* precompute all telemetry-to-mask detectfn values */ 

    int c,m,t,hri;
    double par[3];   
    double r;
    fnptr zfn;
    double hrsum=0.0;
    double normal=1.0;

    zfn = getzfnr(fn);                                          
    for (c=0; c<cc; c++) {
	par[0] = 1; // gsbval[c];
	par[1] = gsbval[cc + c];
	// par[2] = gsbval[2*cc + c];
	if ((fn==14) || (fn==16))
	    //normal[c] = 1/(2 * M_PI * par[0] * par[1] * par[1]);
	    normal = telemscale/(2 * M_PI * par[1] * par[1]);
	else
	    error ("telemetry only coded for HHN and HEX");

	for (m=0; m<mm; m++) {
	    for (t=0; t<nt; t++) {
		r = sqrt(d2 (m, t, mask, xy, mm, nt)); 
		hri = i3(c, m, t, cc, mm); 
		hr[hri] = zfn(par, r) * normal;
                hrsum += hr[hri];   /* test */
	    }	    
	}
    }
    // Rprintf("cc %3d hrmean %10.8f\n", cc, hrsum/mm/nt/cc);
}
/*=============================================================*/

double fnu (int m, int n, int x, int w[], int PIA[], int start[],
	    double hr[], int detect[],
	    int cc, int nc, int kk, int ss, int mm, double minp) {

/* Probability density of telemetry locations for animal n if it belongs to
   latent class x and is centred at m. It is assumed that telemetry data correspond
   to the last detector (kk)

   detspec is an array with the occasion-specific detector type in positions 0:(ss-1)
   and the starting location of fixes for each animal in positions ss:(ss+nc-1).

   The array hr contains precomputed detection function evaluations 
   (dimension cc x mm x nt where nt is total number of telemetry locations).

   The argument 'normal' contains the normalising coefficient 1/c.

   [The size of hr can grow very large, and maybe there needs to be an option to
    compute hr dynamically.]

   Assume detections are sorted by occasion.

 */
    int s,c = 0;
    int hri = 0; 
    int i, t, wi, wxi;
    int count, cumcount = 0;
    double result = 1.0;

    /* fnui <- apply(hr*norm, 1, function (tt) tapply(tt, ID, prod)) */
    /* nsum <- apply(fnui, 1, sum) * (1/mm)   ## area of pixel? */


    for (s=0; s<ss; s++) {
	if (detect[s] == 13) {
	    wi = i3(n, s, kk-1, nc, ss);
	    count = w[wi];  /* number of telemetry fixes */
	    if (count>0) {
		wxi = i4(n, s, kk-1, x, nc, ss, kk);
		c = PIA[wxi] - 1;
		if (c<0) {
		    Rprintf("n %4d m %5d x %5d ss %4d kk %4d wxi %7d PIA[wxi] %4d\n", 
			    n,m,x,ss,kk,wxi,PIA[wxi]);
		    error ("telemetry usage zero on telemetry occasion");
		}
		for (i=cumcount; i<(cumcount+count); i++) {
		    t = start[n] + i;
		    hri  = i3(c, m, t, cc, mm); 
		    result *= hr[hri]; 
		    /* if (result<minp) { */
		    /* 	result = 0.0; */
		    /* 	break; */
		    /* } */
		}
		cumcount += count;
	    }
	}
    }
    return (result);

}
/*=============================================================*/

double prwisignal
   (int m, int n, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
    double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    SIGNAL STRENGTH DETECTOR
    Modified 2012-01-30 to allow missing signal strengths
    Does not need cut; simplified dnorm call 2012-02-01
    Not clear if any benefit in including effort 2012-12-17
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
    double g1;

    spherical = (int) detspec[3];
    for (s=0; s<ss; s++) {
        for (k=0; k<kk; k++) {
            wxi = i4(n, s, k, x, nc, ss, kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                                           /* skip if this trap not set */
                wi = i3(n, s, k, nc, ss);
                count = abs(w[wi]);                                 /* ignore 'deaths */
                if (count==0) {
                    /* not detected at this mic */
                    gi = i3(c,k,m,cc,kk);
		    g1 = gk[gi];
                    result *= countp (0, binomN[s], g1);
                }
                else {
                    /* detected at this mic */
                    start = (int) detspec[wi+4];
                    for (j=0; j<count; j++) {
			sig = signal[start+j];
			if (sig >= 0) {
                            /* valid measurement of signal */
			    /* mu = mufn (k, m, gsbval[c], gsbval[cc + c],
				       traps, mask, kk, mm, spherical); */
			    mu  = mufnL (k, m, gsbval[c], gsbval[cc + c], 
					  dist2, kk, spherical);

			    sdS = gsbval[cc * 2 + c];
			    result *= dnorm((sig - mu), 0, sdS, 0);
			}
			else  {
                            /* signal value missing; detection only */
			    gi = i3(c,k,m,cc,kk);
			    g1 = gk[gi];
			    result *= countp (1, binomN[s], g1);
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
   (int m, int n, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
    double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])

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
    nd = (int) detspec[1];
    spherical = (int) detspec[5];
    /* cut = detspec[2]; */
    muN = detspec[3];
    sdN = detspec[4];
    for (s=0; s<ss; s++) {
        for (k=0; k<kk; k++) {
            wxi = i4(n, s, k, x, nc, ss, kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                                           /* skip if this trap not set */
                wi = i3(n, s, k, nc, ss);
                count = abs(w[wi]);                                 /* ignore 'deaths */
                if (count==0) {
                    /* not detected at this mic */
                    gi = i3(c,k,m,cc,kk);
                    result *= countp (0, binomN[s], gk[gi]);
                }
                else {
                    /* detected at this mic */
                    start = (int) detspec[wi+6];
                    for (j=0; j<count; j++) {
			sig = signal[start+j];
			nois = signal[start+nd+j];

			if (sig >= 0) {
                            /* valid measurement of signal */
			    /* mu = mufn (k, m, gsbval[c], gsbval[cc + c],
				       traps, mask, kk, mm, spherical); */
			    muS  = mufnL (k, m, gsbval[c], gsbval[cc + c], 
					  dist2, kk, spherical);

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
			    result *= countp (1, binomN[s], gk[gi]);
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

double prwipolygon
   (int m, int n, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
    double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])

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
    int nd;
    int start = 0;
    double hint;
    double g1;
    double Tski;

    nd = (int) detspec[1];
    for (s=0; s<ss; s++) {
        if (binomN[s] < 0) error ("negative binomN < 0 not allowed in C fn prwipoly");
        if (pID[s + ss*x]>0)  /* test 2015-11-01 */
        for (k=0; k<kk; k++) {
            wi = i3(n,s,k,nc,ss);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            wxi = i4(n,s,k,x,nc,ss,kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                          /* skip if this polygon not used */
                gi  = i3(c,k,m,cc,kk);
		Tski = Tsk[s * kk + k];
		if (binomN[s]== 0) {                      /* Poisson */
//		    result *= dpois (count, Tski * gk[gi]  * pID[s + ss*x] , 0);   
		    result *= dpois (count, Tski * hk[gi]  * pID[s + ss*x] , 0);   
		}
		else if (binomN[s] == 1) {                 /* Binomial, size from Tsk */ 
		    g1 = gk[gi] * pID[s + ss*x]; 
		    result *= countp (count, round(Tski), g1);
		}
		else {                                  /* Binomial including Bernoulli */ 
		    g1 = gk[gi] * pID[s + ss*x]; 
		    if (fabs(Tski-1) > 1e-10) 
			g1 = 1 - pow(1 - g1, Tski);
		    result *= countp (count, binomN[s], g1);
		}
                /* for each detection pdf(xy) | detected */
                if ((result > minp) && (count>0)) {               /* avoid underflow 2009 11 16 */
		    /* retrieve hint = integral2D(zfn(x) over k)) */
		    hint = hk[gi] / gsbval[c] * detspec[2+c];
                    start = (int) detspec[wi+cc+2];
                    for (j=start; j < start+count; j++) {
                        result *= zfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / hint;
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

/* NEEDS UPDATING FOR 3.0 */
double prwipolygonX
   (int m, int n, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
    double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])

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
    double hint;
    double Tski;

    nd = (int) detspec[1];

    for (s=0; s<ss; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        if (htemp < fuzz) { result = 0; break; }
        wi = nc * s + n;
        k = w[wi];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        /* PR | detected */
        if (k > 0) {
	    Tski = Tsk[s * kk + k-1];
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = Tski * hk[gi];
            pks *= (1-expmin(-htemp)) / htemp;
        }
        /* PR | not detected */
        else  {
            pks = expmin(-htemp);      /* Not captured */
        }
        result *= pks;
        /* Pr(location) */
        if (k > 0) {
	    /* retrieve hint = integral2D(zfn(x) over k)) */
            /* use of  gsbval[c] appropriate for HHN, HHR etc. */
	    hint = hk[gi] / gsbval[c] * detspec[2+c];
            j = (int) detspec[wi+cc+2];
            result *= zfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / hint;
/*
TESTING
	if (m==400 && n == 3)
	    Rprintf("prwipolygonX s %4d htemp %12.8f hk[gi] %12.8f detspec[2+c] %12.8f  hint %12.8f result %12.8f \n", s, 
		    htemp, hk[gi], detspec[2+c], hint, result);
*/
        }

        if (dead) break;
    }
    return (result);
}
/*=============================================================*/

double prwitransect
    (int m, int n, int x, int w[], double xy[], double signal[],
     int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
     int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
     double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])

/*
    Likelihood component due to capture history n (0 <= n < nc)
    given that animal's range centre is at m
    TRANSECT DETECTOR
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of part 0 <= k < kk  */
    int j;   /* index of xy record */
    int c, wxi, wi, gi;
    long count;
    int dead = 0;
    double result = 1.0;
    int nd;
    int start = 0;
    double hint; 
    double g1;
    double Tski;

    nd = (int) detspec[1];
    for (s=0; s<ss; s++) {
        for (k=0; k<kk; k++) {
            wi = i3(n,s,k,nc,ss);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            wxi = i4(n,s,k,x,nc,ss,kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                               /* skip if this transect not used */
		Tski = Tsk[s * kk + k];
                gi  = i3(c,k,m,cc,kk);
		if (binomN[s] == 0) {                      /* Poisson */
//		    result *= countp (count, 0, Tski * gk[gi]); 
		    result *= countp (count, 0, Tski * hk[gi]); 
		}
		else if (binomN[s] == 1) {                 /* Binomial, size from Tsk */ 
		    g1 = gk[gi];
		    result *= countp (count, round(Tski), g1);
		}
		else {                                  /* Binomial, size from binomN */ 
		    g1 = gk[gi];
		    if (fabs(Tski-1) > 1e-10) 
			g1 = 1 - pow(1 - g1, Tski);
		    result *= countp (count, binomN[s], g1);
		}
                if ((result > minp) && (count>0)) {     /* avoid underflow 2009 11 16 */
                    start = (int) detspec[wi+cc+2];
		    /* retrieve integral1D(zfn(x) over k)) */
		    hint = hk[gi] / gsbval[c] * detspec[2+c];
                    for (j=start; j < start+count; j++) {
                        result *= zfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / hint;
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

double prwitransectX
   (int m, int n, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], double hk[], int binomN[], double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, zfnptr zfn, double gsbval[],
    double traps[], double dist2[], double Tsk[], double mask[], double minp, double pID[])
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
    double hint;
    double Tski;

    nd = (int) detspec[1];

    for (s=0; s<ss; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        if (htemp < fuzz) { result = 0; break; }
        wi = nc * s + n;
        k = w[wi];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        /* PR | detected */
        if (k > 0) {
	    Tski = Tsk[s * kk + k-1];
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = Tski * hk[gi];
            pks *= (1-expmin(-htemp)) / htemp;
        }
        /* PR | not detected */
        else  {
            pks = expmin(-htemp);      /* Not captured */
        }
        result *= pks;
        /* Pr(location) */
        if (k > 0) {
            j = (int) detspec[wi+cc+2];
            hint = gk[gi] / gsbval[c] * detspec[2+c];
            result *= zfn(j, m, c, gsbval, cc, xy, mask, nd, mm, 0) / hint;
        }

        if (dead) break;
    }
    return (result);
}
/*=============================================================*/

void pdotpoint (double *xy, int *nxy, double *traps, double *dist2, 
		int *detect, double *Tsk, int *markocc, int *ss, int *kk, int *fn, double *par, 
		double *w2, int *binomN, double *value)
{
    int i,k,s;
    double dk2;
    double tempval;
    double g0;
    double sigma;
    double z = 1;
    double cutval [1];
    double p;
    double Tski = 1.0;
    cutval[0] = 0;
    int allsighting = 1;
    
    for (s=0; s<*ss; s++) {
        if (markocc[s]>0) allsighting = 0;                    /* capture occasions */
    }
    /*-----------------------------------------------------------*/
    if (dist2[0] < 0) {
        if (anypolygon(detect, *ss) || anytransect(detect, *ss)) {
	    dist2 = (double *) S_alloc(1, sizeof(double));
	}
	else {
	    dist2 = (double *) S_alloc(*kk * *nxy, sizeof(double));
	    makedist2 (*kk, *nxy, traps, xy, dist2);
	}
    }
    else {
	squaredist(*kk, *nxy, dist2);
    }
    /*-----------------------------------------------------------*/

    if (*fn>18)
        error("pdotpoint requires detectfn < 18");
    g0 = par[0];
    sigma = par[1];
    if (!((*fn == 0) || (*fn == 2) || (*fn == 4) || (*fn == 9) 
	  || (*fn == 14) || (*fn == 16)))
        z = par[2];
    if ((*fn == 10) || (*fn == 11) || (*fn == 12) || (*fn == 13))
        cutval[0] = par[3];
    
    for (i=0; i<*nxy; i++) {
        tempval = 1;
        for (s=0; s<*ss; s++) {
            if (((markocc[s]>0) || allsighting) && (detect[s]!=13)) {                     
                for (k=0; k<*kk; k++) {
                    Tski = Tsk[s * *kk + k];
                    if (Tski > 1e-10) {
                        dk2 = d2L(k, i, dist2, *kk);
                        p = pfn(*fn, dk2, g0, sigma, z, cutval, *w2);
                        /* counts */
                        if (detect[s] == 2) {    
                            if (binomN[s] == 0)
                                /* 2015-12-22  */
                                p = 1 - countp(0, 0, Tski * hazard(p));
                            else if (binomN[s] == 1)
                                p = 1 - countp(0, round(Tski), p);
                            else {
                                if (fabs(Tski-1) > 1e-10)
                                    p = 1 - pow(1-p, Tski);
                                p = 1 - countp(0, binomN[s], p);
                            }
                        }
                        else {                 /* multi, proximity */ 
                                if (fabs(Tski-1) > 1e-10)
                                    p = 1 - pow(1-p, Tski);
                        }
                        tempval *= 1 - p;
                    }
                }
            }
        }
        value[i] = 1 - tempval;
    }
}

/*=============================================================*/

// 2017-04-04 NOT TESTED YET

void hdotpoint (double *xy, int *nxy, double *traps, double *dist2, 
		int *detect, double *Tsk, int *occasions, int *ss, int *kk,
                int *fn, double *par, 
		double *w2, double *value)
{
    int i,k,s;
    double dk2;
    double lambda0;
    double sigma;
    double z = 1;
    double h, H = 0;
    double Tski = 1.0;
    double miscparm[1];
    
    /*-----------------------------------------------------------*/
    if (dist2[0] < 0) {
        if (anypolygon(detect, *ss) || anytransect(detect, *ss)) {
	    dist2 = (double *) S_alloc(1, sizeof(double));
	}
	else {
	    dist2 = (double *) S_alloc(*kk * *nxy, sizeof(double));
	    makedist2 (*kk, *nxy, traps, xy, dist2);
	}
    }
    else {
	squaredist(*kk, *nxy, dist2);
    }
    /*-----------------------------------------------------------*/
    if ((*fn<14) || (*fn>18))
        error("hdotpoint requires detectfn 14-18");
    lambda0 = par[0];
    sigma = par[1];
    if (!((*fn == 14) || (*fn == 16))) z = par[2];
    for (i=0; i<*nxy; i++) {
        H = 0.0;
        for (s=0; s<*ss; s++) {
            if (occasions[s]>0) {                     
		for (k=0; k<*kk; k++) {
		    Tski = Tsk[s * *kk + k];
		    if (Tski > 1e-10) {
			dk2 = d2L(k, i, dist2, *kk);
			h = zfn(*fn, dk2, lambda0, sigma, z, miscparm, *w2);
			H += Tski * h;
		    }
		}
	    }
        }
        value[i] = H;
    }
}

/*=============================================================*/

void hdotpoly (double *xy, int *nxy, double *traps, int *detect, double *Tsk, int *markocc, 
               int *nk, int *ss, int *kk, int *fn, double *par, double *value)
{
    int i,k,s;
    double hk = 0;
    double sumhk;
    int *cumk;
    double H = 1.0;
    double *ex;
    double Tski = 1.0;
    ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
    int allsighting = 1;
    
    for (s=0; s<*ss; s++) {
        if (markocc[s]>0) allsighting = 0;                    /* capture occasions */
    }
    
    cumk = (int *) R_alloc(*nk+1, sizeof(int));
    cumk[0] = 0;
    for (i=0; i<*nk; i++)
        cumk[i+1] = cumk[i] + kk[i];

    /* 2017-03-22 integrated hazard */
    if (anypolygon(detect, *ss))           /* polygon */
	H = hintegral(*fn, par);     
    else if (anytransect(detect, *ss))     /* transect */
	H = hintegral1(*fn, par);
    else 
	error ("unrecognised detector type in hdotpoly");
    
    for (i=0; i<*nxy; i++) {
        sumhk = 0.0;
        for (s=0; s<*ss; s++) {
            if (((markocc[s]>0) || allsighting) && (detect[s]!=13)) {                     
                for (k=0; k<*nk; k++) {
                    Tski = Tsk[s * *nk + k];
                    if (Tski > 1e-10) {  /* effort > 0 */
                        if (anypolygon(detect,*ss))
                            hk = par[0] * integral2D (*fn, i, 0, par, 1, traps, xy,
                                               cumk[k], cumk[k+1]-1, cumk[*nk], *nxy, ex) / H;
                        else if (anytransect(detect,*ss))
                            hk = par[0] * integral1D (*fn, i, 0, par, 1, traps, xy,
                                               cumk[k], cumk[k+1]-1, cumk[*nk], *nxy, ex) / H;
                        sumhk += hk * Tski;
                    }
                }
            }
        }
        value[i] = sumhk;
    }
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

void getdenom (int *fn, double *miscparm, double *mask, int *mm, double *scale,
	       double sigma, double z) {
    /* fill third column of mask with 1 / denom[m] or zero if no usage*/
    /* mask must be dim(*mm, 3) if miscparm[1]==0, otherwise dim(*mm,4) */
    /* accepts interrupts once every checkinterval mask rows */

    double sigma2, lam, dsq;
    double d = 1.0;
    double expzj = 1.0;
    double denomtot = 0.0;
    int m, j, im;
    int checkinterval = 100;
    /*---------------------------------------------------------*/
    /* normalization 2013-11-10 */
    if (((*fn == 4) || ((*fn >= 14) && (*fn <= 18))) && (fabs(miscparm[0]) > 0.5)) {
	sigma2 =  sigma * sigma;
	for (m=0; m<*mm; m++) {
	    im = 2 * *mm + m;
	    mask[im] = 0;
	    for (j=0; j<*mm; j++) {
		if ((fabs(miscparm[1]) > 0.5))
		    expzj = mask[3 * *mm + j];  /* exp(covariate val) at mask point j */
		dsq = d2(j, m, mask, mask, *mm, *mm);
		if (*fn != 14) d = sqrt(dsq);
		if (*fn == 14) lam = exp(-dsq / 2 / sigma2);
		else if (*fn == 4) lam = (double) (d <= sigma);
		else if (*fn == 15) lam = 1 - exp(- pow(d /sigma , - z));
		else if (*fn == 16) lam = exp(-d / sigma);
		else if (*fn == 17) lam = exp(-(d-z)*(d-z) / 2 / sigma2);
		else if (*fn == 18) lam = pgamma(d,z,sigma/z,0,0);
		else error("unrecognised fn");
		mask[im] += lam * expzj;
	    }
	    denomtot += mask[im];
	    if (m % checkinterval == 0)
		R_CheckUserInterrupt();
	}
        /* optional scaling */
	*scale = 1;
	if (fabs(miscparm[2]) > 0.5) {
	    *scale = denomtot / *mm;
	    for (m = 0; m < *mm; m++) {
		im = 2 * *mm + m;
		if (mask[im] > 0)
		    mask[im] = *scale / mask[im];
	    }
	}
    }	
}
/*=============================================================*/

void getdenomext (int *fn, double *miscparm, double *mask, int *mm, double *sigma, 
		  double *z, double *invdenom, double *scale) {
    /* external wrapper for getdenom : 1 / denom[m] or zero if no usage*/
    /* mask must be dim(*mm, 3) if miscparm[1]==0, otherwise dim(*mm,4) */
    /* accepts interrupts once every checkinterval mask rows */
    int m, im;
    getdenom (fn, miscparm, mask, mm, scale, *sigma, *z);
    for (m = 0; m < *mm; m++) {
	im = 2 * *mm + m;
	invdenom[m] = mask[im];
    }	
}
/*=============================================================*/

void precompute(
    int detect[], 
    int fn, 
    int binomN[],
    int ss,
    int kk,
    int mm,
    int cc,
    int nk,
    int cumk[],
    double traps[],
    double dist2[],   /* added 2014-08-27 */
    double mask[],
    double gsbval[],
    double miscparm[],
    double detspec[],
    double gk[],
    double hk[],
    int debug
 ) {
    /*---------------------------------------------------------*/
    /* populate pre-computed gk and gk0, hk and hk0  arrays    */
    /*
        detect[s] may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  exclusive polygon detector
        4  exclusive transect detector
        5  signal detector
        6  polygon detector
        7  transect detector
     	12 signalnoise
	13 telemetry
    */
    /*---------------------------------------------------------*/

    int c,m,k,gi;

    double par[4];   /* passing parameter values to integr fn  */
    double H = 1;
    gfnLptr gfnL;
    double *ex; 
    gfnL = getgfnL(fn); 

    /* allowtelem = 1 */
    if (allpoint(detect, ss, 1, 1)) {
	for (c=0; c<cc; c++) {
	    /*---------------------------------------------------------*/
	    /* normalization NOT implemented 2013-11-10 */
	    /* getdenom(fn, miscparm, mask, mm, scale, gsbval[*cc + c], gsbval[2 * *cc + c]); */
	    /*---------------------------------------------------------*/
	    for (k=0; k<nk; k++) {
		for (m=0; m<mm; m++) {
                    gi = i3(c,k,m,cc,nk);
                    gk[gi] = gfnL(k, m, c, gsbval, cc, dist2, 
				  nk, miscparm);
                    hk[gi] = -log(1-gk[gi]);   /* stopgap */
                }
            }
        }
    }
    else if (anypolygon(detect, ss)) {
	ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
	for (c=0; c<cc; c++) {
	    /*---------------------------------------------------------*/
	    /* normalization NOT implemented 2013-11-11 */
	    /* getdenom(fn, miscparm, mask, mm, scale, gsbval[*cc + c], gsbval[2 * *cc + c]); */
	    /*---------------------------------------------------------*/
	    par[0] = gsbval[c];
	    par[1] = gsbval[cc + c];
	    par[2] = gsbval[2 * cc + c];
	    H = hintegral(fn, par);         /* 2017-03-22 unbounded integrated hazard */

	    if (debug>2) Rprintf("c %4d  H %12.8f\n", c, H);

	    detspec[2+c] = H;               /* passed to prwipolygon */
	    for (k=0; k<nk; k++) {               /* over parts */
		for (m=0; m<mm; m++) {
		    gi = i3(c, k, m, cc, nk);
		    /* 2017-03-22 strictly use hazard form : expected detections of animals at m */
                    /* par[0] only makes sense here if fn is HHN, HHR, HEX, HAN HCG */
		    hk[gi] = par[0] * integral2D (fn, m, 0, par,
						  1, traps, mask, cumk[k], cumk[k+1]-1, cumk[nk],
						  mm, ex) / H;
		    gk[gi] = 1 - exp(-hk[gi]);
/*TESTING
 if (m==400) Rprintf("hk[gi]  %12.8f\n", hk[gi]);
*/		}
	    }
	}
    }
    else if (anytransect(detect,ss)) {
	ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
	for (c=0; c<cc; c++) {
	    /*---------------------------------------------------------*/
	    /* normalization NOT implemented 2013-11-11 */
	    /* getdenom(fn, miscparm, mask, mm, scale, gsbval[*cc + c], gsbval[2 * *cc + c]); */
	    /*---------------------------------------------------------*/
	    par[0] = gsbval[c];
	    par[1] = gsbval[cc + c];
	    par[2] = gsbval[2* cc + c];
	    H = hintegral1(fn, par);         /* 2017-03-22 unbounded integrated hazard */
	    detspec[2+c] = H;               /* passed to prwitransect */
	    for (k=0; k<nk; k++) {               /* over transects */
		for (m=0; m<mm; m++) {
		    gi = i3(c,k,m,cc,nk);
		    /* 2017-03-22 strictly use hazard form */
		    hk[gi] = par[0] * integral1D (fn, m, c, gsbval, cc,
			 traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], mm, ex) / H;
		    gk[gi] = 1 - exp(-hk[gi]);
		}
	    }
	}
    }	
    else {
	Rprintf("allpoint %12d \n", allpoint(detect, ss, 1, 1));
        error ("unrecognised detector type in external C fn");
    }

} /* end precompute */
/*=============================================================*/

prwfnptr getprwfn (int detect[], int ss) {
/* ss unused? */
    prwfnptr prwfn;
    prwfn = prwipoint; /* default for detect 0,1,2,13 */
//    if (allmulti(detect,ss))
//	prwfn = prwimulti; /* old for 0 */
//    else 
    if (detect[0] == 3)
        prwfn = prwipolygonX;
    else if (detect[0] == 4)
        prwfn = prwitransectX;
    else if (detect[0] == 5)
	prwfn = prwisignal;
    else if (detect[0] == 6)
        prwfn = prwipolygon;
    else if (detect[0] == 7)
        prwfn = prwitransect;
    else if (detect[0] == 12) {
    	prwfn = prwisignalnoise;   /* experimental 2012-02-07 */
    }
    return(prwfn);
} /* end getprwfn */
/*=============================================================*/

int getstart(int detect[], int start[], int nc1, int nc, 
		 int ss, int nk, int w[]) {

	/*------------------------------------------------------------*/
	/* Identify start positions of ancillary data for each animal */
        /* Telemetry is incompatible with polygon or signal detectors */

	int nd = 0;
	int i,s,k,wi,first;

	if (anytelemetry(detect, ss)) {
            /* assume all telemetry fixes associate with last detector */
	    for (i=0; i< nc; i++) {
		first = 1;
		for (s=0; s<ss; s++) {
		    if (first) start[i] = nd;
		    first = 0;
		    wi = i3(i,s,nk-1,nc,ss);
		    nd += abs(w[wi]);
		}
	    }
	}
	else {
	    if (anypolygon(detect, ss) || anytransect(detect,ss) || 
		anysignal(detect,ss)) {
		/* start[z] indexes the first row in xy (or element in signal)
		   for each possible count z, where z is w-order (isk) */
		for (k=0; k<nk; k++) {
		    for (s=0; s< ss; s++) {
			for (i=0; i< nc; i++) {
			    wi = i3(i,s,k,nc,ss);
			    start[wi] = nd;
			    nd += abs(w[wi]);
			}
		    }
		}
	    }
	}
	return(nd);
    }
/*=============================================================*/

void getdetspec (int detect[], int fn, int nc,  int nc1, 
		 int cc, int nmix, int nd, int nk, int ss, int kk, int mm, 
		 int PIA[], double miscparm[], int start[], double detspec[]) {

    /* detector-type-specific data passed later to prwi functions */

    /* mixtures are group-specific for full likelihood, and     */
    /* individual-specific for conditional likelihood           */

    int i,s;

    /* default for (possibly mixed) point detectors and telemetry */
    if (allpoint(detect, ss, 0, 1)) {
	for(s=0; s<ss; s++)
	    detspec[s] = (double) detect[s];

	/* start position (in xy) of telemetry fixes of each animal */
	if (anytelemetry(detect, ss)) {
	    for (i=0; i< nc; i++)
		detspec[ss+i] = (double) start[i];
	}
    }
    /*  polygonX and transectX */
    else if ((detect[0] == 3) || (detect[0] == 4)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        /* pass index of first detection of each animal + occasion */
        /* maximum 1 per polygon as exclusive */
        for (i=0; i< (nc * ss); i++)
            detspec[2+cc+i] = (double) start[i];
    }
    /*  polygon and transect */
    else if ((detect[0] == 6) || (detect[0] == 7)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        /* pass index of first detection of each animal + occasion */
        /* maximum 1 per polygon as exclusive */
        for (i=0; i< (nc * ss * nk); i++)
            detspec[2+cc+i] = (double) start[i];
    }
    /* signal */ 
    else if (detect[0] == 5) {    
        for (i=0; i<3; i++) detspec[i]= miscparm[i];
        detspec[3]= ((fn == 11) || (fn == 13));     /* spherical */
        for (i=0; i< (nc * ss * nk); i++)
            detspec[4+i] = (double) start[i];
    }
    /* signal-noise */
    else if (detect[0] == 12) {                           
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;                     /* number of detections or detectors??*/
	detspec[2]= miscparm[0];                      /* cut */
	detspec[3]= miscparm[1];                      /* noise mean */
	detspec[4]= miscparm[2];                      /* noise sd */
        detspec[5]= ((fn == 11) || (fn == 13));     /* spherical */
        for (i=0; i< (nc * ss * nk); i++)
            detspec[6+i] = (double) start[i];
    }
}
/*=============================================================*/

void geth (int nc1, int cc, int nmix, int nk, int ss, int mm, 
	   int PIA[], int hc0[], double hk[], double Tsk[],
		 double h[], int hindex[]) {

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

    int fullns = 0;
    int next = 0;
    int c,c0,i,m,n,s,k,x,gi,hi;
    int PIAval0, PIAvalk;   /* added 2012-12-12 */
    double Tski;            /* added 2012-12-17 */

    /* Recognise when not fully specified by n.s, c  */
    /* this arises when model = bk, Bk               */
    /* and when detector covariates vary by time     */
    /* Modified for speed 2012-12-12 by excluding    */
    /* case that detectors are not used PIA < 0      */
    
    for (n=0; n < nc1; n++) {
        for (s=0; s < ss; s++) {
	    PIAval0 = PIA[i4(n,s,0,0, nc1, ss, nk)];
	    for (k=1; k<nk; k++) {
                PIAvalk = PIA[i4(n,s,k,0, nc1, ss, nk)];
		if (PIAval0 < 0)
		    PIAval0 = PIAvalk;
		else if (PIAvalk>0) {
		    if (PIAval0 != PIAvalk) {
			fullns = 1;
			break;
		    } 
		}              
	    } 
	    if (fullns) break;
	}
	if (fullns) break;
    }
    /* Rprintf("fullns %5d \n", fullns); */

    for (i=0; i<cc; i++) hc0[i] = -1;
    next = 0;        
    for (n=0; n < nc1; n++) {
	for (s=0; s < ss; s++) {
	    hi = s*nc1 + n;
	    /* Case 1. within-trap variation */
	    if (fullns) {
		for (k=0; k < nk; k++) {
		    Tski = Tsk[s * nk + k];
		    for (x = 0; x < nmix; x++) {
			c = PIA[i4(n,s,k,x, nc1, ss, nk)]-1; 
			for (m = 0; m < mm; m++) { 
			    if (c >= 0) {
				gi = i3(c,k,m,cc,nk);
				h[i3(x,m,hi,nmix, mm)] += Tski * hk[gi];
			    }
			}
		    }
		}
		hindex[hi] = hi;   
	    }
	    /* Case 2. no within-trap variation */
	    else {
		c0 = PIA[i4(n,s,0,0, nc1, ss, nk)] - 1;                    
		if (hc0[c0] < 0) {
		    hc0[c0] = next;
		    next ++;
		    for (k=0; k < nk; k++) {
			Tski = Tsk[s * nk + k];
			for (x = 0; x < nmix; x++) {
			    c = PIA[i4(n,s,k,x, nc1, ss, nk)]-1; 
			    for (m=0; m< mm; m++) { 
				if (c >= 0) {
				    gi = i3(c,k,m,cc,nk);
				    h[i3(x,m, hc0[c0],nmix, mm)] += Tski * hk[gi];
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
/*=============================================================*/

/* no need to return gpar 2015-10-10 */

void getpmix(int gpar, int nc1, int nmix, int knownclass[], int nc, int cc, int ss, 
            int nk, int grp[], int PIA[], double gsbval[], double pmixg[], 
            double pmixn[]){

    /*-------------------------------*/
    /* mixture proportions           */
    /* by group and by animal        */

    int g, n, c, x, wxi;
    double pmix;
    if (nmix>1) {
        /* one extra real parameter for h2, h3 */
        gpar++;
        for (n=0; n<nc1; n++) {
            for (x=0; x<nmix; x++) {
                wxi = i4(n,0,0,x,nc,ss,nk);
                c = PIA[wxi] - 1;
		if (c<0) error ("c<0 error in getpmix");  /* 2017-02-08 */
		pmix = gsbval[cc * (gpar-1) + c];  /* last column in gsbval */

                /* group-specific, and overall pmix by class for knownclass case */
		g = grp[n]-1;
		pmixg[nmix * g + x] = pmix;

                /* individual-specific */                
		if (knownclass[n] > 1) {
		    if (knownclass[n] == (x+2))   /* knownclass=2 maps to x=0 */
			pmixn[nmix * n + x] = 1;
		    else 
			pmixn[nmix * n + x] = 0;
		}
		else
		    pmixn[nmix * n + x] = pmix;
            }
        }
    }
}
/*=============================================================*/

int nval(int detect0, int nc1, int cc, int ss, int nk) {
    /* compute and allocate the space needed for ancillary data */
    int nval;
    if ((detect0==3) || (detect0==4))
        nval = 2 + cc + nc1 * ss;
    else if ((detect0==6) || (detect0==7))
        nval = 2 + cc + nc1 * ss * nk;
    else if (detect0==5)
        nval = 4 + nc1 * ss * nk;
    else if (detect0==12)    /* signalnoise */
        nval = 6 + nc1 * ss * nk;
    else
        nval = ss + nc1;    /* default for (possibly mixed) 0,1,2,13 */
    return(nval);
}
/*=============================================================*/

int fillcumk(int detect[], int ss, int kk[], int cumk[]){
    /* determine number of polygons if polygon detector */
    /* for polygon detectors, kk is vector ending in zero */
    /* and first nk+1 elements of vector cumk are filled */
    int i;
    int nk = 0;
    if (anypolygon(detect, ss) || anytransect(detect, ss)) {
        cumk[0] = 0;
        for (i=0; i<maxnpoly; i++) {
            if (kk[i]<=0) break;
            cumk[i+1] = cumk[i] + kk[i];
            nk++;
        }
	for (i=0; i<nk; i++)                /* over parts */
	    if ((cumk[i+1] - cumk[i]) > maxvertices)
		error("exceeded maximum number of vertices %d per polygon", maxvertices);
    }
    else
        nk = kk[0];  /* return number of detectors unchanged */

    return(nk);
}
/*=============================================================*/

void fillng(int nc, int gg, int grp[], int ng[]) {
    int g, n;
    /* Count number per group (not used for CL)                */
    /* Assume histories sorted by group = individual           */
    /* CH are numbered 0 <= n < nc in C code                  */
    for (g=0; g<gg; g++)
        ng[g] = 0;
    for (n=0; n<nc; n++) { 
        g = grp[n] - 1; 
        ng[g]++;
    }
}
/*=============================================================*/

int filla0 (int like, int n, int markocc[], int ncol, int PIA0[],
	    double gk0[], double hk0[], int detect[], int binomN[], double Tsk[], int ss, int nk, 
	    int mm, int cc0, int nmix, double gsb0val[], int allsighting, 
	    double pimask[], double area, double a0[]) {
    double pdt;
    int x;
    int m;
    if (pimask[0] < 0) return(5);

    for (x=0; x < nmix; x++) {
	if (like == 6) {
	    a0[x] = 0;
	    for (m=0; m< mm; m++) {
		pdt = pndot (m, 0, markocc, x, ncol, PIA0, gk0, hk0, detect, binomN, 
			     Tsk, ss, nk, cc0, nmix, gsb0val, allsighting);
		/* adjust both integral and pi to size of pixel */
		a0[x] += pimask[m] * pdt * area * mm;  
	    }
	    if (a0[x] < 0) {    // changed from <= 2015-12-31
		// Rprintf("bad a0 in filla0\n");
		return(6);
	    }
	}
	else a0[x] = 1; 
    }
    /* resultcode 5,6 if bad */
    return(0);   // successful completion
}
/*==============================================================================*/

void integralprw1 (
    int    *detect,    /* detector 0 multi, 1 proximity etc. */
    double *gsb0val,   /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *grp,       /* group number for 0<=n<*nc   [full likelihood only] */
    int    *nc,        /* number of individuals */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    int    *gg,        /* number of groups */
    int    *nmix,      /* number of mixtures */
    int    *knownclass,  /* known membership of 'latent' classes */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *dist2,     /* distances (optional: -1 if unused) */
    double *Tsk,       /* nk x s usage matrix */
    int    *markocc,   /* which are marking occasions? */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *cc0,       /* number of g0/sigma/b combinations [naive animal] */
    int    *PIA0,      /* lookup which g0/sigma/b combination to use for given n, S, K
                          [naive animal] */
    int    *ncol,      /* number of columns in PIA0; added 2009 06 25 */
    double *area,      /* area associated with each mask point (ha) */
    double *miscparm,  /* miscellaneous parameters (cutval, normalization, NT etc) */
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
    int i,n,m,s,x;
    double asum = 0;
    double *gk0 = NULL;
    double *hk0 = NULL;
    double D = 1.0;
    int cumk[maxnpoly];
    int nk = 0;
    int nc1;
    double *pmixg;      /* proportion in each mixture by group*/
    double *pmixn;      /* proportion in each mixture by individual*/
    int gpar;
    double *detspec = NULL;
    int nv;
    double *pID = NULL;
    int allsighting = 1;

    *resultcode = 1;                   /* generic failure */

    nk = fillcumk(detect, *ss, kk, cumk);                            /* detections per polygon */

    pID = (double *)  S_alloc(*ss * *nmix, sizeof (double));
    gpar = par3(*fn) + 2;   /* number of detection parameters */
    gpar = markresightini (*ss, *nmix, markocc, nk, *ncol, PIA0, 
			   *cc0, gsb0val, pID, gpar);
    // compute allsighting rather than base on like >= 5
    for (s=0; s< *ss; s++)
	if (markocc[s] > 0) allsighting = 0;

    /*---------------------------------------------------------*/
    /* Adjust for zero animals                                 */
    /* (otherwise fails with nc = 0 because *nc also used to   */
    /* dimension arrays etc. so replace nc1 = max(1, *nc)      */
    if (*nc == 0) 
        nc1 = 1;
    else
        nc1 = *nc;
    /*---------------------------------------------------------*/

    nv = nval(detect[0], nc1, *cc0, *ss, nk);

    detspec = (double *) R_alloc(nv, sizeof(double));

    /* Allocate space for array of naive detection probability */
    gk0 = (double *) R_alloc(*cc0 * nk * *mm, sizeof (double));
    hk0 = (double *) R_alloc(*cc0 * nk * *mm, sizeof (double));

    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    pmixn = (double *)  R_alloc(*nc * *nmix, sizeof (double));
    for (i=0; i < *gg * *nmix; i++) pmixg[i] = 1; /* default */
    for (i=0; i < nc1 * *nmix; i++) pmixn[i] = 1; /* default */
    getpmix(gpar, nc1, *nmix, knownclass, *nc, *cc0, *ss, nk, 
		   grp, PIA0, gsb0val, pmixg, pmixn);

    if (dist2[0] < 0) {
	dist2 = (double *) S_alloc(*kk * *mm, sizeof(double));
	makedist2 (*kk, *mm, traps, mask, dist2);
    }
    else {
	squaredist(*kk, *mm, dist2);
    }

    precompute(detect, *fn, binomN, *ss, *kk, *mm, *cc0, nk, cumk,   
	       traps, dist2, mask, gsb0val, miscparm, detspec, gk0, hk0, 0); 


    for (n=0; n < nc1; n++) {            /* CH numbered 0 <= n < nc1 */
        if ((*ncol > 1) || (n == 0)) {   /* no need to repeat if constant */
            if ((n+1) > *ncol) {         /* groups not implemented */
                *resultcode = 3;
                return;
            }
            asum = 0;
            for (x=0; x<*nmix; x++) {
                for (m=0; m<*mm; m++) {
                    if (*useD > 0) D = mask[2 * *mm + m];
                    asum += pmixn[*nmix * n + x] * D * 
			pndot (m, n, markocc, x, *ncol, PIA0, gk0, hk0, detect, 
			       binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, 
			       allsighting);
                }
            }
        }
        a[n] = *area * asum;
    }
    *resultcode = 0;                   /* successful completion */
}

/*==============================================================================*/

/*  check if we need to consider variation among individuals */
/* i.e. check if detection parameters constant for given s,k */
int anyvarying (
    int    nc,     /* number of capture histories (or groups if PIA0 has that dim) */
    int    ss,     /* number of occasions */
    int    nk,     /* number of traps */
    int    nmix,   /* number of mixture classes */
    int    *PIA0    /* lookup which g0/sigma/b combination to use for given n, S, K [naive] */
) {
    int i,n,s,k,x;
    int wxi;
    int indiv = 0;
    for (s=0; s<ss; s++) {
	for (k=0; k<nk; k++) {
	    for (x=0; x<nmix; x++) {
		wxi = i4(0,s,k,x,nc,ss,nk);        /* changed from *kk 2011-02-07 */
		i = PIA0[wxi];
		for (n=1; n<nc; n++) {
		    wxi = i4(n,s,k,x,nc,ss,nk);    /* changed from *kk 2011-02-07 */
		    if (i != PIA0[wxi]) {
			indiv = 1; break;
		    }
		}
	    }
	}
    }
    return(indiv);
}
/*==============================================================================*/

clock_t timestamp(clock_t ticks1, int *counter) {
    clock_t ticks2;
    ticks2 = clock();
    Rprintf("check %2d: time used %12d ticks\n", *counter, ticks2-ticks1); 
    *counter = *counter + 1;
    return(ticks2);
}
/*==============================================================================*/

double classmembership (int n, int x, int knownclass[], double pmixn[], int nmix) {

/* Return probability individual n belongs to class x. This may be binary 
   (0/1) in the case of known class, or continuous if class is unknown */

    double pmixnx = 0;
    int knownx = -1;

    if (knownclass[n] == 1) 
	knownx = -1;                         /* unknown class */
    else
	knownx = imax2(0, knownclass[n]-2);  /* known class */
    
    /* unknown class : weighted by probability of membership  */
    /* known class does not require probability of membership */
    if (knownx < 0)
	pmixnx = pmixn[nmix * n + x];
    else if (knownx == x)
	pmixnx = 1.0;
    else 
	pmixnx = 0.0;
    return pmixnx;

}

void zerok (int detect[], int nc, int ss, int nk, int w[]) {
    /* Bring any capture to k=0 for exclusive-detector occasions */
    /* It saves time to do this once, rather than for each mask point */
    /* even better to do this once before calling C */
    int n, s, k, w0;
    int wi = 0;

    for (n=0; n<nc; n++) {
	for (s=0; s<ss; s++) {
	    if (detect[s]==0 || detect[s]==3 || detect[s]==4) {
		k = 0;
		do {
		    wi = i3(n, s, k, nc, ss);
		    k++;
		}
		while ((w[wi] == 0) && (k<nk));
		    
		if (w[wi] != 0) {                              /* Captured */
		    w0 = i3(n, s, 0, nc, ss);
		    if (w[wi]<0) w[w0] = -k;           /* save detector number */
		    else w[w0] = k; 
		    // Rprintf("n %5d s %5d  w[w0] %5d  \n", n, s, w[w0]);
		}
	    }
	}
    }
}
/*=============================================================*/

void sumwsk (int w[], int nc, int ss, int nk, double musk[]) {
    int n,s,k;
    for (n=0;n<nc; n++)
	for (s=0; s<ss; s++)
	    for(k=0; k<nk; k++)
		musk[k + s*nk] += abs(w[i3(n,s,k,nc,ss)]); 
}
/*=============================================================*/

void getnmarked(int w[], int nc, int ss, int nk, int nmarked[]) {  // Number marked at occasion s
    int n,s,k;
    int *caught;
    caught = (int *)  S_alloc(nc, sizeof (int)); 
    for (s=0; s<ss; s++) {
	for (n=0;n<nc; n++) {
	    if (!caught[n]) {
		for(k=0; k<nk; k++) {
		    if (abs(w[i3(n,s,k,nc,ss)]) > 0) {
			nmarked[s] ++; 
			caught[n] = 1;
			break;
		    }
		    if (caught[n]) break;
		}
	    }
	}
	// Rprintf("s %4d nmarked %4d \n", s, nmarked[s]);
	nmarked[s+1] = nmarked[s];  // requires ss+1 slots
    }
}
/*=============================================================*/

void selectoccasions (int ss, int markocc[], int detect[], double detspec[], int debug, int mcode) {
    int s;
    /* reset if mcode < -5 */
    if (debug>=2) Rprintf("detspec ");
    for(s=0; s<ss; s++) {
	if ((markocc[s]==mcode) || (mcode < -5))
	    detspec[s] = (double)detect[s];
	else 
	    detspec[s] = -1.0;
	if (debug>=2) Rprintf(" %4.0f", detspec[s]);
    }
    if (debug>=2) Rprintf("\n");
}
/*=============================================================*/

void secrloglik (
    int    *like,        /* likelihood 0 full, 1 conditional etc. See list at top */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
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
    int    *knownclass,  /* known membership of 'latent' classes; 1='unknown' */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *dist2,       /* distances (optional: -1 if unused) */
    double *Tsk,         /* nk x s usage matrix */
    int    *telemcode,   /* 0 none 1 IT 2 DT 3 CT */
    int    *markocc,     /* which are marking occasions? */
    int    *Tu,          /* detector x occasion matrices of no. sightings Tu */
    int    *Tm,          /* detector x occasion matrices of no. sightings Tm */
    int    *Tn,          /* detector x occasion matrices of no. sightings Tn */
    double *chat,        /* sighting overdispersion */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *Dmask,       /* density at each point on mask, possibly x group */
    double *pimask,      /* individual probability density; used if pimask[0] >= -tol (>=0) */
    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    double *gsb0val,     /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *cc0,         /* number of g0/sigma/b combinations for naive animals */
    int    *PIA,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *PIA0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive] */
    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameters (cutval, normalization, NT etc.) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *minprob,     /* minimum value of P(detection history) */
    double *telemscale,
    int    *debug,       /* 0 messages off, 1 messages on */
    double *comp,        /* likelihood components */
    double *value,       /* returned log likelihood (sum of components) */
    int    *resultcode   /* 0 if OK */
)

/*
  A note on ordering of input data arrays 2012-02-09

  Using 'i' to subscript animals, 's' to subscript occasions, 'k' to subscript detectors..

 'w' is in isk order - dim(CH) = c(nc,ss,nk) 

 'signal' is in linear ksi order and includes positive detections only, 1:nd 

 For signal-noise detectors (*detect==12) the positions nd+1 : 2nd are noise measurements 

 The integer array start (dim nc x ss x nk or nc x ss) holds the index for each isk to
   the first corresponding detection in 'signal' 

 It is possible in principle for there to be more than one detection per isk; these will
   follow in sequence, hence the name 'start' 

*/

{
    /* indices */
    int    i,n,g,m,x;

    /* miscellaneous */
    double temp, tempsum, tempp, templog, prwi;
    size_t memrequest;   // "size_t is defined in stddef.h which the header defining R_alloc..."

    /* group arrays */
    int    *ng;       /* number per group */
    int    *nm;       /* number per known mixture class */
    double *sumD;     
    double *sumDp;    
    double *pdot = NULL;

    /* generalised detection and detector functions */
    zfnptr  zfn;
    prwfnptr prwfn;

    /* numbers of individuals and detectors */
    int    nc1;
    int    cumk[maxnpoly];
    unsigned long nk = 0;

    /* number of 'g' (detection) parameters */
    int    ncol;
    int    gpar;    

    /* number of detections */
    int    nd = 0;  

    /* number of marking occasions and related */
    double *pID = NULL;
    double *pdots = NULL;
    double Tulik = 0.0;
    double Tmlik = 0.0;
    double Tnlik = 0.0;
    int    allsighting = 0;
    double numerator = 1;
    double a0[maxnmix];
    double *Tumu = NULL;
    double *Tmmu = NULL;
    double *tmpmusk = NULL;

    int    *firstocc = NULL;
    int    alltelem = 0;
    int    anytelem = 0;
    int    *ntelem;    /* number telemetry histories per group */
    int    *telem;     /* telemetry history? */
    int    *nzero;     /* number all-zero histories per group  */
    int    *kcnzero;   /* number all-zero histories per knownclass  */
    int    *allzero;   /* all-zero history? */
    int     nonzero;   

    /* switch between Tm likelihoods: 0 shorter (unconditional) 1 longer (conditional) */
    int Tmswitch = 1;      
    int *nmarked = NULL;      

    /* stored pdot for allsighting */
//    double *pdot = NULL;

    /* pre-computed detection probabilities */
    double *gk = NULL;
    double *gk0 = NULL;
    double *hk = NULL;
    double *hk0 = NULL;
    double *hr = NULL;

    /* mixture membership probability */
    double *pmixg = NULL;
    double *pmixn = NULL;
    double pmixnx;        /* proportion in current class */
    int    known = 0;     /* number of 'known-class' individuals */
    int    knownx = -1;   /* class of current animal (negative if unknown) */

    /* passing data to prwfn */
    double *detspec = NULL;
    int    *start = NULL;

    /* total hazard computation */
    int    *hc0 = NULL;        
    int    *hindex = NULL;     
    double *h = NULL;          

    /* CL only */
    int    indiv = 0;     /* indicate detection varying between individuals */
    double asum[maxnmix] = {0.0};  /* initialises all to zero */
    double *a = NULL;

    double pdt = 0;
    int    nv;
    double tempN = 0;
    double p;
    double tol = 1e-8;    /* numerical tolerance */

    /* telemetry pi */
    double f = 1.0;         /* probability density of range centre */

    /*-------------------------------------------------------*/

    /* MAINLINE */

     clock_t ticks;
     int timing = 0;
     int counter = 0;
     ticks = clock();
     for (i=0; i<6; i++) comp[i]=0;
     *value = 0;

     /* 0 ******************************************************/
     if (timing) ticks = timestamp(ticks, &counter);

    /*---------------------------------------------------------*/
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
    /*---------------------------------------------------------*/

    /* Number per known mixture, and total in known class */
    /* On input, knownclass 1 is used for 'unknown' : */
    /* knownclass 1 -> nm[0] 'unknown' */
    /* knownclass 2 -> nm[1] 'latent class 1' */
    /* knownclass 3 -> nm[2] 'latent class 2' */
    nm = (int *) R_alloc(*nmix + 1, sizeof(int));
    kcnzero = (int *) S_alloc(*nmix+1, sizeof(int));
    fillng(*nc, *nmix+1, knownclass, nm);                           
    if (*nmix > 1) {
	for (x = 1; x < (*nmix+1); x++)
	known += nm[x];
    }
    /*---------------------------------------------------------*/

    /* Number per group */
    ng = (int *) S_alloc(*gg, sizeof(int));
    ntelem = (int *) S_alloc(*gg, sizeof(int));
    nzero = (int *) S_alloc(*gg, sizeof(int));
    telem = (int *) S_alloc(*nc, sizeof(int));
    allzero = (int *) S_alloc(*nc, sizeof(int));
    fillng(*nc, *gg, grp, ng);                                    
    pdot = (double *)  S_alloc(nc1 * *nmix * *mm, sizeof (double));   // limited use: CT?

    /*---------------------------------------------------------*/

    /* Detections per polygon cumk and nk */
    nk = fillcumk(detect, *ss, kk, cumk);                            
    /*---------------------------------------------------------*/

    /* Select detection function and detector */
    /* see utils.c for these functions */
    /* zfn is hazard function used by polygon and transect detectors */
    zfn = getzfn(*fn);                                          
    prwfn = getprwfn(detect, *ss);
    gpar = par3(*fn) + 2;

    /*---------------------------------------------------------*/

    /* Mark-resight */
    /* Update pID, gpar */
    pID = (double *)  S_alloc(*ss * *nmix, sizeof (double));
    if (*like == 1) ncol = nc1;
    else ncol = *gg;
    gpar = markresightini (*ss, *nmix, markocc, nk, ncol, PIA, *cc, gsbval,
			   pID, gpar);
    allsighting = (*like >= 5);
    if (Tu[0]>=0) {
        /* storage for Pr(marked on or before s) */
	pdots = (double *)  S_alloc(*ss * *nmix * *mm, sizeof (double));
    }
    Tumu  = (double *)  S_alloc(*ss * nk, sizeof (double));  // Tu expected value
    Tmmu = (double *)  S_alloc(*ss * nk, sizeof (double));  // Tm expectedvalue
    nmarked = (int *)  S_alloc(*ss+1, sizeof (int));  // Number marked at occasion s
    if (Tm[0]>=0 && Tmswitch) {
	tmpmusk = (double *)  S_alloc(*ss * nk, sizeof (double));
    }

    if (chat[0]<=0 || chat[1]<=0 || chat[2]<=0) {
	Rprintf("chat must be positive\n");
	*resultcode = 11;
	return;
    }

    // Number marked at occasion s; for Tu distrib = binomial
    if ((Tu[0]>=0) && *distrib) {
        getnmarked(w, *nc, *ss, nk, nmarked);  
    }

    /*---------------------------------------------------------*/
    alltelem = alltelemetry(detect,*ss);
    anytelem = anytelemetry(detect,*ss);

    /*---------------------------------------------------------*/
    /* 'start' is index to position of ancillary data for each
       individual, used for polygon, transect and sound detectors */
    if ( (detect[0] == 3) || (detect[0] == 4) ) 
	start = (int *) R_alloc(nc1 * *ss, sizeof(int));
    else if (((detect[0]>=5) && (detect[0]<=9)) || (detect[0]==12)) 
	start = (int *) R_alloc(nc1 * *ss * nk, sizeof(int));
    else if (anytelem)
	start = (int *) R_alloc(nc1, sizeof(int));
    nd = getstart(detect, start, nc1, *nc, *ss, nk, w);
    if (*debug>2) Rprintf("nd %4d \n", nd);

    /* work arrays for telemetry */
    if (anytelem) {
	memrequest =  *cc;
	memrequest *= *mm;
	memrequest *= nd;
	if (*debug>=2) {
	    Rprintf("Memory request: cc %6d mm %6d nd %6d cc.mm.nd %lu \n",
		    *cc, *mm, nd, memrequest);
	}
	hr = (double *) S_alloc(memrequest, sizeof (double));
    }

    /*---------------------------------------------------------*/

    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    pmixn = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    for (i=0; i < (*gg * *nmix); i++) pmixg[i] = 1; /* default */
    for (i=0; i < (nc1 * *nmix); i++) pmixn[i] = 1; /* default */
    getpmix(gpar, nc1, *nmix, knownclass, *nc, *cc, *ss, nk, 
		   grp, PIA, gsbval, pmixg, pmixn);
    if (*debug>1) Rprintf("pmixg %6.4f, %6.4f\n", pmixg[0], pmixg[1]);

    nv = nval(detect[0], nc1, *cc, *ss, nk);
    detspec = (double *) R_alloc(nv, sizeof(double));

    memrequest =  *cc;
    memrequest *= *mm;
    memrequest *= nk;
    if (*debug>=2) {
	Rprintf("Memory request: cc %6d nk %6d mm %6d cc.nk.mm %lu \n",
		*cc, nk, *mm, memrequest);
    }
    gk = (double *) S_alloc(memrequest, sizeof(double)); 
    hk = (double *) S_alloc(memrequest, sizeof(double)); 
    memrequest = *cc0;
    memrequest *= *mm;
    memrequest *= nk;
    gk0 = (double *) S_alloc(memrequest, sizeof(double));
    hk0 = (double *) S_alloc(memrequest, sizeof(double));

    /* 1 *********************************************************/
    if (timing) ticks = timestamp(ticks, &counter);

    /*-----------------------------------------------------------*/
    if (dist2[0] < 0) {
        if (anypolygon(detect,*ss) || anytransect(detect,*ss)) {
	    dist2 = (double *) S_alloc(1, sizeof(double));
	}
	else {
	    dist2 = (double *) S_alloc(*kk * *mm, sizeof(double));
	    makedist2 (*kk, *mm, traps, mask, dist2);
	}
    }
    else {
	squaredist(*kk, *mm, dist2);
    }
    // dist2 checked 2017-01-01
    /*----------------------------------------------------------*/

    precompute(detect, *fn, binomN, *ss, *kk, *mm, *cc, nk, cumk,    
	       traps, dist2, mask, gsbval, miscparm, detspec, gk, hk, *debug); 
    if (anyb(gsbval, gsb0val, *cc, gpar)) {
	precompute(detect, *fn, binomN, *ss, *kk, *mm, *cc0, nk, cumk,   
		   traps, dist2, mask, gsb0val, miscparm, detspec, gk0, hk0, *debug); 
    }
   else {   // merely copy 
	for (i=0; i<(*cc * *mm * nk); i++) {
	    gk0[i] = gk[i];
	    hk0[i] = hk[i];
	}
        // detspec[2+c] remains unchanged 
    }

    firstocc = (int *) R_alloc (nc1, sizeof(int));
    if (anytelem)  {
	gethr(detect, *fn, *ss, *mm, *cc, nd, mask, xy, gsbval, hr, *telemscale);
    }
    /* getfirstocc2(*ss, nk, *nc, w, grp, detect, firstocc, ntelem, telem, nzero, allzero); */
    getfirstocc2(*ss, nk, *nc, w, grp, knownclass, detect, firstocc, ntelem, telem, 
		 nzero, kcnzero, allzero);
    if (*debug>2) for (x=0;x<(*nmix+1);x++) Rprintf("x %4d kcnzero[x] %5d\n", x, kcnzero[x]);
    if (*debug>2) Rprintf("ntelem[0] %4d \n", ntelem[0]);

    /* 2 *********************************************************/
    if (timing) ticks = timestamp(ticks, &counter);
    /*-----------------------------------------------------------*/

    R_CheckUserInterrupt();

    if (anyexclusive(detect,*ss)) {   /* exclusive detectors require hazard */ 
	zerok (detect, *nc, *ss, nk, w);
        hc0 = (int *) R_alloc (*cc, sizeof(int));
        hindex = (int *) S_alloc (nc1 * *ss, sizeof(int));
	h = (double *) S_alloc (nc1 * *ss * *mm * *nmix, sizeof(double)); 
	geth (nc1, *cc, *nmix, nk, *ss, *mm, PIA, hc0, hk, Tsk, h, hindex);
    } 

   /* complete filling of detspec */
    getdetspec (detect, *fn, *nc, nc1, *cc, *nmix, nd, nk, *ss, 
		*kk, *mm, PIA, miscparm, start, detspec);

    if (*debug>2)  {
	Rprintf("Detector ");
	for(i=0;i<*ss;i++) Rprintf("%3d", (int) detspec[i]);
	Rprintf("\n");
    }

    /*-------------------------*/
    /* Now evaluate likelihood */
    /*-------------------------*/

    if (*like == 1) {  /* Conditional likelihood */
        indiv = anyvarying (*nc, *ss, nk, *nmix, PIA0);
	if (*debug>2) Rprintf("indiv %3d\n", indiv);
	if (!indiv && !alltelem) {
	    /* all individuals the same */
	    /* save time by doing this once, rather than inside n loop */
	    for (x=0; x<*nmix; x++) {
		asum[x] = 0;
		for (m=0; m<*mm; m++) {
		    asum[x] += pndot (m, 0, markocc, x, *nc, PIA0, gk0, hk0, detect, binomN, 
				 Tsk, *ss, nk, *cc0, *nmix, gsb0val, allsighting);

		    /* for mark-resight likelihood */
		    if (Tm[0] >= 0) {
			error("unfinished");
		    }
		}
	    }
	}
	/* else asum not needed or calculated for each individual in loop below */
	if (*debug>1) Rprintf("asum(0) %10.6f\n", asum[0]);

	/* 3 *********************************************************/
	if (timing) ticks = timestamp(ticks, &counter);
        /*-----------------------------------------------------------*/

	/* Loop over individuals... */
	a = (double *)  S_alloc(nc1, sizeof (double));
	for (n=0; n<*nc; n++) {                      /* CH numbered 0 <= n < *nc */
	    a[n] = 0;
	    if (knownclass[n] == 1) 
		knownx = -1;                         /* unknown class */
	    else
		knownx = imax2(0, knownclass[n]-2);  /* known class */

	    /* zero array for expected sightings of unidentified */
	    if (Tm[0] >= 0 && Tmswitch) {
		for(i=0; i< (*ss * nk); i++) tmpmusk[i] = 0;
	    }

	    /* unknown class : weighted by probability of membership  */
	    /* known class does not require probability of membership */
	    tempsum = 0;

	    for (x = 0; x < *nmix; x++) {
		temp = 0;
		if (knownx < 0)
		    pmixnx = pmixn[*nmix * n + x];
		else if (knownx == x)
		    pmixnx = 1.0;
		else 
		    pmixnx = 0.0;

		if (pmixnx > 1e-6) {
		    if (indiv > 0)
			asum[x] = 0;
		    for (m=0; m<*mm; m++) {
			if (indiv && *telemcode<3 && !allzero[n]) { 
			    pdt = pndot (m, n, markocc, x, *nc, PIA0, gk0, hk0, detect, 
					 binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, 
					 allsighting); 
			    asum[x] += pdt;
			}

			/* telemetry fixes IT, DT, CT */
                        /* could provide non-mask-based likelihood for IT */
			if (*telemcode>0 && telem[n]) {
			    f = fnu (m, n, x, w, PIA, start, hr, detect,  
				      *cc, *nc, *kk, *ss, *mm, *minprob);
			}
			else {
			    f = 1.0;
			}

			if (f > *minprob) 
			{   /* default f = 1 ; minprob was tol*/
			    if (telem[n] && *telemcode == 1)
				prwi = 1.0;
			    else
				prwi = prwfn (m, n, x, w, xy, signal, PIA, gk, hk, binomN,
			    		  detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix,
			    		  zfn, gsbval, traps, dist2, Tsk, mask, *minprob, pID);
			    temp += prwi * f;
			    if (*debug>1 &&  m==400) 
				Rprintf("m %5d n %4d f %12.10f prwi %12.10f sum prwi %12.10f\n",
					m,n,f,prwi,temp);
			    if (Tm[0]>=0  && Tmswitch) { 
				incmusk (temp, n, m, x, w, PIA, gk, hk, binomN, *cc, *nc, nk, *ss,
					 *nmix, gsbval, Tsk, markocc, firstocc[n], detect, tmpmusk);
			    }
			}
		    }
		    a[n] += asum[x] * pmixnx;  // weight by pmixnx 2017-01-02
		    tempsum += temp * pmixnx;
		    if (*debug>2) Rprintf("%5d tempsum %12.8f a[n] %9.6f \n", n, tempsum, a[n]);
		}
	    }    /* end of loop over mixture classes */

            if (Tm[0]>=0 && Tmswitch) finmusk (*ss, nk, tmpmusk, Tmmu, tempsum, pID, 1); 

            /* 2017-02-07 let templog go NaN if it wants... */
	    templog = log(tempsum * *area);
	    if (!R_FINITE(templog)) *resultcode = 9;
	    if (*resultcode == 9) return;
	    comp[0] += templog;

            /* if all telemetry then no adjustment for p(detected) */
	    if (!allzero[n] && (*telemcode<3 || !telem[n])) 
		comp[1] -= log(a[n] * *area);

	    // Rprintf("n %5d value %12.8f\n", n, comp[0]);
	    R_CheckUserInterrupt();
	}        /* end of loop over individuals */

        /* multinomial probability of known class membership (excludes coefficient) */
        /* ignored if telemetry */
	if (known>0) {
	    tempsum = 0;
	    for (x=0; x<*nmix; x++) {
		asum[x] = 0;
		for (m=0; m<*mm; m++) {
		    asum[x] += pndot (m, 0, markocc, x, *nc, PIA0, gk0, hk0, detect, 
				      binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, 
				      allsighting);
		}
		tempsum += asum[x] * pmixg[x];
	    }
	    for (x=0; x<*nmix; x++) 
		comp[2] +=  (nm[x+1] - kcnzero[x+1]) * log(asum[x] * pmixg[x] / tempsum); 
	    if (!R_FINITE(comp[2])) *resultcode = 9;
    	    if (*resultcode == 9) return;
	}

	/* sightings of unidentified animals */
	if (Tm[0]>=0) {
	    if (!Tmswitch) {
                /* standalone calculation of Tmmu */
		*resultcode = expectedTmTu (*like, *distrib, 0, *nc, *ss, nk,
					  *cc0, *nmix, pmixg, *mm, Dmask, pimask,
					    *area, markocc, pdots, ncol, PIA0, gk0, hk0, 
					    binomN, detect, Tsk, nmarked, asum, pID, Tmmu);
		if (*resultcode>0) return;
	    }
	    *resultcode = Tsightinglik (Tm, *ss, nk, markocc, ncol, detect, Tsk, Tmmu,
					 *debug, &Tmlik);
	    if (*debug>=1) {
		Rprintf("CL sighting marked nonID Tmlik %12.8f \n", Tmlik);  
	    }
	    if (*resultcode>0) return;
	    comp[5] = Tmlik / chat[1];
	}
	
    } /* end (*like == 1) */
    /*-------------------------------------------------------------------------------------------*/

    else if ((*like == 0) || (*like == 2)) {  /* Full likelihood; 2 = split sightings */

	sumD = (double *) S_alloc(*gg, sizeof(double));
	sumDp = (double *) S_alloc(*gg, sizeof(double));
	if ((known>0) && (*nmix>1)) {
	    /* no groups, only mixture classes; use group 0 for aggregates */
	    for (m=0; m<*mm; m++)  {
		sumD[0] += Dmask[m];
		for (x=0; x<*nmix; x++) {
		    pdt = pndot (m, 0, markocc, x, *gg, PIA0, gk0, hk0, detect, binomN,
				 Tsk, *ss, nk, *cc0, *nmix, gsb0val, allsighting);
		    sumDp[0] += pdt * pmixg[x] * Dmask[m]; 
		}
	    }
	}
	else {
	    /* groups; incompatible with mixtures containing some known */
	    for (g=0; g<*gg; g++) {
		for (m=0; m<*mm; m++)  {
		    sumD[g] += Dmask[*mm * g + m];
		    for (x=0; x<*nmix; x++) {
			pdt = pndot (m, g, markocc, x, *gg, PIA0, gk0, hk0, detect, binomN, 
				     Tsk, *ss, nk, *cc0, *nmix, gsb0val, allsighting);
			sumDp[g] += pdt * pmixg[*nmix * g + x] * Dmask[*mm * g + m];
			pdot[i3(m,g,x, *mm, *gg)] = pdt; // save if needed for CT?
			/* save time- and site-specific Pr(marked) */
                        /* for mark-resight likelihood */
			if (Tu[0] >= 0) {
			    getpdots (m, g, markocc, x, *gg, PIA0, gk0, hk0, detect,
				    binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, pdots);
			}
		    }
		}
	    }
	}
	
	/* 3 *********************************************************/
	if (timing) ticks = timestamp(ticks, &counter);
        /*-----------------------------------------------------------*/

        
        /*------------------------------------------*/
        /* compute likelihood component from pr(wi) */
        /*------------------------------------------*/
	g = 0;

	for(i=0; i< *ss; i++) 
	    if (markocc[i] == -1) detspec[i] = -2.0;


	for (n=0; n < *nc; n++) {
	    tempsum = 0;
	    if ((known==0) || (*nmix==1)) 
		g = grp[n]-1;
	    if (Tm[0]>=0  && Tmswitch) {
		for(i=0; i< (*ss * nk); i++) tmpmusk[i] = 0;
	    }
	    if (*like == 2) {
		selectoccasions(*ss, markocc, detect, detspec, *debug, 1);
	    }
	    for (x=0; x<*nmix; x++) {
		pmixnx = classmembership (n, x, knownclass, pmixn, *nmix);
		if (pmixnx > 1e-6) {
		    for (m=0; m<*mm; m++) {
			/* telemetry fixes */
			if (*telemcode > 0 && telem[n]) {
			    f = fnu (m, n, x, w, PIA, start, hr, detect,
				      *cc, *nc, *kk, *ss, *mm, *minprob);
			}
			else f = 1.0;

			if (f > *minprob) {
			    if (telem[n] && *telemcode == 1)  /* independent telem */
				prwi = 1.0;
			    else
				prwi = prwfn (m, n, x, w, xy, signal, PIA, gk, hk, binomN,
			    		  detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix,
			    		  zfn, gsbval, traps, dist2, Tsk, mask, *minprob, pID);
			    tempp = f * prwi * pmixnx;
                            /* do not do this for IT, allzero CT telemetry animals */	    
			    if (!telem[n] || !allzero[n])
				tempp *= Dmask[*mm * g + m];  /* Dprwi; note pdt cancels except CT (see next) */	 
			    if (*telemcode==3 && !allzero[n] && telem[n]) {
				tempp *= pdot[i3(m,g,x, *mm, *gg)];
			    }
   
			    tempsum += tempp;                             /* sumDprwi */
			    /* for unidentified marked sightings, increment each mu_sk */
			    if (Tm[0]>=0 && Tmswitch) { 
				incmusk (tempp, n, m, x, w, PIA, gk, hk, binomN, *cc, *nc, nk, *ss,
					 *nmix, gsbval, Tsk, markocc, firstocc[n], detect, tmpmusk);
			    }
			}
		    }
		}
	    }
            /* for unidentified marked sightings */	        
	    if (Tm[0]>=0 && Tmswitch) {
		finmusk (*ss, nk, tmpmusk, Tmmu, tempsum, pID, 1);		
	    }
	    templog = log(tempsum);
	    if (!R_FINITE(templog)) *resultcode = 9;
	    if (*resultcode == 9) return;
	    if (*debug>=2) Rprintf("log prwi1: n %4d logprwi %8.3f \n", n, templog);
	    comp[0] += templog;

	    /*--------------------------------------------------------*/
	    /* separate treatment of sighting component */
	    /* later be more careful here - maybe should include in pmix loop above */
	    if (*like == 2) {
		selectoccasions(*ss, markocc, detect, detspec, *debug, 0);
		tempsum = 0;
		for (x=0; x<*nmix; x++) {
		    pmixnx = classmembership (n, x, knownclass, pmixn, *nmix);
		    if (pmixnx > 1e-6) {
			for (m=0; m<*mm; m++) {
			    if (f > tol) {    /* default f = 1 */
				prwi = prwfn (m, n, x, w, xy, signal, PIA, gk, hk, binomN, 
					      detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, 
					      zfn, gsbval, traps, dist2, Tsk, mask, *minprob, pID);
				tempp = f * prwi * pmixnx;	    
				tempsum += tempp;                                /* sumprwi */
			    }
			}
		    }
		}	        
		templog = log(tempsum);
		if (!R_FINITE(templog)) *resultcode = 9;
		if (*resultcode == 9) return;
		if (*debug>=2) Rprintf("log prwi2: n %4d logprwi %8.3f \n", n, templog);
		selectoccasions(*ss, markocc, detect, detspec, *debug, -10);  /* reset */
		comp[0] += templog;
	    }
	    /*--------------------------------------------------------*/

	    R_CheckUserInterrupt();
	} /* end n loop over indiviudals */

	/* 4 *************************************************************/
	if (timing) ticks = timestamp(ticks, &counter);	    
	/*---------------------------------------------------------------*/

	/* Known membership of latent classes                            */
	/*---------------------------------------------------------------*/

	if ((known>0) && (*nmix>1)) {
	    nonzero = *nc - nzero[0];         

	    /*-------------------------------*/
	    /* likelihood component due to n */
	    /*-------------------------------*/
	    /* Poisson */
	    if (*distrib == 0) {
		comp[1] = gpois(nonzero, sumDp[0] * *area, 1);
	    }
	    /* binomial */
	    if (*distrib == 1) {
		tempN = sumD[0] * *area;
		if (nonzero > tempN) {
		    *resultcode = 8;
		    return;
		}
		p = sumDp[0] / sumD[0];
		comp[1] = dbinom_raw (nonzero, tempN, p, 1-p, 1);
	    }
	    /*---------------------------------------------------*/
	    /* likelihood component due to denominator of f()    */
	    comp[2] -= nonzero * log(sumDp[0]);

	    /* adjustment for mixture probabilities when class known */
	    for (x=0; x<*nmix; x++) {
		comp[3] += (nm[x+1] - kcnzero[x+1]) * log(pmixg[x]);
	    }
	}

	/*---------------------------------------------------------------*/
	/* No known membership of latent classes                         */
	/*---------------------------------------------------------------*/
	else {
	    for (g=0; g<*gg; g++) {

        	if (*telemcode==0 || *telemcode==2)  /* no telem or DT */
		    nonzero = ng[g];
		else                                 /* IT or CT */
		    nonzero = ng[g] - nzero[g];         

		/*-------------------------------*/
		/* likelihood component due to n */
		/*-------------------------------*/
	    
		/* Poisson */
		if (*distrib == 0) {
			comp[1] += gpois(nonzero, sumDp[g] * *area, 1);
		}
		/* binomial */
		else if (*distrib == 1) {
		    tempN = sumD[g] * *area;
		    if (nonzero > tempN) {    /* including ntelem[g] */
			*resultcode = 8;
			return;
		    }
		    p = sumDp[g] / sumD[g];
		    comp[1] += dbinom_raw (nonzero, tempN, p, 1-p, 1);
		}
		
		/*---------------------------------------------------*/
		/* likelihood component due to denominator of f()    */
		/*---------------------------------------------------*/
		comp[2] -= nonzero * log(sumDp[g]); 
	    }
	}

	/*------------------------------------------------------*/
	/* sightings of unmarked or unidentified animals        */
	/* array pdots previously filled by function getpdots   */
	/*------------------------------------------------------*/

	if ((Tu[0]>=0) || (Tm[0]>=0))
	    if (*gg > 1) error("mark-resight does not allow groups");

	/*------------------------------------------------------*/
	if (Tu[0]>=0) {	    /* unmarked */
	    *resultcode = expectedTmTu (*like, *distrib, 1, *nc, *ss, nk,
					*cc0, *nmix, pmixg, *mm, Dmask, pimask, *area, 
					markocc, pdots, ncol, PIA0, gk0, hk0, binomN, 
					detect, Tsk, nmarked, asum, pID, Tumu);
	    if (*resultcode>0) return;
	    *resultcode = Tsightinglik (Tu, *ss, nk, markocc, ncol, detect, Tsk, Tumu,
					*debug, &Tulik);
	    if (*resultcode>0) return;
	    comp[4] = Tulik / chat[0];
	}
	/*------------------------------------------------------*/
	if (Tn[0]>=0) {	   /* unresolved */
            /* using Tumu for Tnmu  */
            /* using chat[0] for Tn */
	    *resultcode = expectedTmTu (*like, *distrib, 2, *nc, *ss, nk,
					*cc0, *nmix, pmixg, *mm, Dmask, pimask, *area, 
					markocc, pdots, ncol, PIA0, gk0, hk0, binomN, 
					detect, Tsk, nmarked, asum, pID, Tumu);
	    if (*resultcode>0) return;
	    *resultcode = Tsightinglik (Tn, *ss, nk, markocc, ncol, detect, Tsk, Tumu,
					*debug, &Tnlik);
	    if (*resultcode>0) return;
	    comp[4] += Tnlik / chat[2];  /* sum of any Tu, Tn likelihood components */
	}
	/*------------------------------------------------------*/
	if (Tm[0]>=0) {	   /* marked, not identified */
	    if (!Tmswitch) {
		*resultcode = expectedTmTu (*like, *distrib, 0, *nc, *ss, nk,
					    *cc0, *nmix, pmixg, *mm, Dmask, pimask, *area, 
					    markocc, pdots, ncol, PIA0, gk0, hk0, binomN, 
					    detect, Tsk, nmarked, asum, pID, Tmmu);
		if (*resultcode>0) return;
	    }
	    *resultcode = Tsightinglik (Tm, *ss, nk, markocc, ncol, detect, 
					Tsk, Tmmu, *debug, &Tmlik);
	    if (*resultcode>0) return;
	    comp[5] = Tmlik / chat[1];
	}
	/*------------------------------------------------------*/

	if (*debug>=1) {
	    Rprintf("Full likelihood: Total %8.3f,  Tulik %8.3f  chat(Tu) %7.3f",
		    comp[4], Tulik, chat[0]); 
	    Rprintf(" Tmlik %8.3f  chat(Tm) %8.3f \n", Tmlik, chat[1]);
	    Rprintf(" Tnlik %8.3f  chat(Tn) %8.3f \n", Tnlik, chat[2]);
	}

    }   /* end (*like == 0) || (*like == 2) */
    /*----------------------------------------------------------------------------------------*/

    else if ((*like == 5) || (*like == 6)) /* All-sighting conditional likelihood  2015-10-15 */
    {
        /* Preliminaries ... */
        /* all-sighting data presumes prior marking with known distribution over the habitat mask */
/*
	if (*like == 6) {
	    for (x=0; x<*nmix; x++) {
		asum[x] = 0;
		for (m=0; m<*mm; m++) {
		    asum[x] += pndot (m, 0, markocc, x, *nc, PIA0, gk0, hk0, detect, binomN, 
				      Tsk, *ss, nk, *cc0, *nmix, gsb0val, allsighting);
		}
	    }
	}
*/
        /* 3 *********************************************************/
	if (timing) ticks = timestamp(ticks, &counter);
        /*-----------------------------------------------------------*/

	*resultcode = filla0 (*like, 0, markocc, ncol, PIA0, gk0, hk0, detect, binomN, Tsk, 
			      *ss, nk, *mm, *cc0, *nmix, gsb0val, allsighting, pimask, 
			      *area, a0);
	if (*resultcode>0) return;   // abort this likelihood evaluation with resultcode 5 or 6

	if (Tm[0] >= 0 && Tmswitch) {
	    for(i=0; i< (*ss * nk); i++) tmpmusk[i] = 0;
	}
	if (*debug>1)
	    Rprintf("nc = %4d\n", *nc);
	/* Loop over individuals... */
	for (n=0; n<*nc; n++) {                      /* CH numbered 0 <= n < *nc */
	    tempsum = 0;
	    if (Tm[0] >= 0 && Tmswitch) {
		for(i=0; i<(*ss * nk); i++) tmpmusk[i] = 0;
	    }
	    for (x = 0; x < *nmix; x++) {
		temp = 0;
		pmixnx = classmembership (n, x, knownclass, pmixn, *nmix);
		/* Rprintf("pmixnx %8.6f\n", pmixnx); */
		if (pmixnx > 1e-6) {
		    for (m=0; m< *mm; m++) {
			numerator = pmixnx * pimask[m] * 
			    prwfn (m, n, x, w, xy, signal, PIA, gk, hk, binomN, 
				    detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix,
				    zfn, gsbval, traps, dist2, Tsk, mask, *minprob, pID);
			if ((Tm[0]>=0) && Tmswitch && (firstocc[n] < *ss)) { 
			    incmusk (numerator, n, m, x, w, PIA, gk, hk, binomN, *cc, *nc, nk, *ss,
					 *nmix, gsbval, Tsk, markocc, firstocc[n], detect, tmpmusk);
			}
			if (*like == 5)
			    temp += numerator;
			else if (a0[x]>0)
			    temp += numerator / a0[x];  
		    }
		}
		tempsum += temp;
	        /* Rprintf("%5d tempsum %12.8f asum[x] %9.6f \n",
		   n, tempsum, asum[x]);   */
	    }    /* end of loop over mixture classes */

/* MOVE INSIDE X LOOP TO ALLOW FOR NMIX>1 */
	    if (Tm[0] >= 0 && Tmswitch) {
		finmusk (*ss, nk, tmpmusk, Tmmu, tempsum, pID, 1);		
	    }

	    templog = log(tempsum * *area);    
	    if (!R_FINITE(templog)) *resultcode = 9;
	    if (*resultcode == 9) return;
	    comp[0] += templog;
	    R_CheckUserInterrupt();
	}        /* end of loop over individuals */

        /* separately compute for nzero all-zero detection histories */
        /* of pre-marked animals. There are nzero of these observed when */
        /* like == 5; the number is estimated when like == 6 */

	if ((Tm[0]>0) && Tmswitch && (*like == 6)) { 
	    error("allsighting Tm not ready yet");
/*
	    tempsum = 0.0;   
            for(i=0; i<(*ss * nk); i++) tmpmusk[i] = 0;
	    for (x = 0; x < *nmix; x++) {
		temp = 0;
		pmixnx = 1 / *nmix;    // TEMPORARY FUDGE - NEED PMIX 
		if (pmixnx > 1e-6) {
		    for (m=0; m< *mm; m++) {
			numerator = pmixnx * pimask[m] * prwi0(m,x,PIA etc.);
			temp += numerator / a0[x];  
			incmusk (numerator, n, m, x, w, PIA, gk, hk, binomN, *cc, *nc, nk, *ss,
					 *nmix, gsbval, Tsk, markocc, 0, detect, tmpmusk);
		    }		    
		}
		tempsum += temp;
	    } 
	    if (*like == 5)
		finmusk (*ss, nk, tmpmusk, Tmmu, tempsum, pID, nzero);		
	    else 
		finmusk (*ss, nk, tmpmusk, Tmmu, tempsum, pID, *nc/a0[0] - *nc);
	    // note preceding fudge when we don't know a0[x] 
	    */
	}

    /* sighting matrices */

	if (Tu[0]>=0) {	
	    *resultcode = expectedTmTu (*like, *distrib, 1, *nc, *ss, nk,
					*cc0, *nmix, pmixg, *mm, Dmask, pimask, *area, 
					markocc, pdots, ncol, PIA0, gk0, hk0, binomN, 
					detect, Tsk, nmarked, a0, pID, Tumu);
	    if (*resultcode>0) return;
	    *resultcode = Tsightinglik (Tu, *ss, nk, markocc, ncol, detect, Tsk, Tumu,
					*debug, &Tulik);

	    if (*resultcode>0) return;
	    comp[4] = Tulik / chat[0];
	}
	if (Tm[0]>=0) {	
	    if (!Tmswitch) {
	        *resultcode = expectedTmTu (*like, *distrib, 0, *nc, *ss, nk,
                                     *cc0, *nmix, pmixg, *mm, Dmask, pimask, *area, 
                                     markocc, pdots, ncol, PIA0, gk0, hk0, binomN, 
                                     detect, Tsk, nmarked, a0, pID, Tmmu);
	        if (*resultcode>0) return;
	    }
	    *resultcode = Tsightinglik (Tm, *ss, nk, markocc, ncol, detect, Tsk, Tmmu,
                                 *debug, &Tmlik);
                                 if (*resultcode>0) return;
                                 comp[5] = Tmlik / chat[1];
	}
	/* if (Tn[0] > 0) error ("unresolved sightings incompatible with all-sighting"); */
	/* NOTE: Tn not an option for sighting only */
        /* BUT we may use it for one session among several 2017-04-17 */
	if (Tn[0]>=0) {	   /* unresolved */
	    /* using Tumu for Tnmu  */
	    *resultcode = expectedTmTu (*like, *distrib, 2, *nc, *ss, nk,
					*cc0, *nmix, pmixg, *mm, Dmask, pimask, *area, 
					markocc, pdots, ncol, PIA0, gk0, hk0, binomN, 
					detect, Tsk, nmarked, asum, pID, Tumu);
	    if (*resultcode>0) return;
	    *resultcode = Tsightinglik (Tn, *ss, nk, markocc, ncol, detect, Tsk, Tumu,
					*debug, &Tnlik);
	    if (*resultcode>0) return;
	    comp[4] += Tnlik / chat[2];  /* sum of any Tu, Tn likelihood components */
	}
	if (*debug>=1) {
	    Rprintf("Full likelihood: Total %8.3f,  Tulik %8.3f  chat(Tu) %7.3f",
		    comp[4], Tulik, chat[0]); 
	    Rprintf(" Tmlik %8.3f  chat(Tm) %8.3f ", Tmlik, chat[1]);
	    Rprintf(" Tnlik %8.3f  chat(Tn) %8.3f \n", Tnlik, chat[2]);
	}
    } /* end (*like == 5) || (*like == 6) */
    /*----------------------------------------------------------------------------------------*/

    /* General wrap-up - no changes to likelihood */
    for (i=0; i<6; i++) {
	*value += comp[i];  /* sum components */
    }
    if (!R_FINITE(*value)) {
	*resultcode = 9;
	*value = -1e10;
    }
    else {
	*resultcode = 0;   /* successful termination secrloglik */
    }
}

/*===========================================================================================*/
/*===========================================================================================*/
    
/*
    'naive' functions are used to estimate auto initial values
    these use only the halfnormal detection function 4/5/08
*/

void naived (
  double *sigma,   /* Parameter : detection scale */
  int    *kk,      /* number of traps */
  int    *nc,
  int    *wt,      /* integer trap weights */
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
    error ("invalid detection function in external function naived");

    for (n=0; n<*nc; n++)
    {
        x = animals[n];
        y = animals[n + *nc];

        for (i=0; i<*kk; i++)
	    if (wt[i] > 0) {
		for (j=0; j<(i-1); j++) {
		    dij = (traps[i] - traps[j]) * (traps[i] - traps[j]) +
                        (traps[i+*kk] - traps[j+*kk]) * (traps[i+*kk] - traps[j+*kk]);
		    d21 = (traps[i] - x) * (traps[i] - x) + (traps[i+*kk] - y) * (traps[i+*kk] - y);
		    d22 = (traps[j] - x) * (traps[j] - x) + (traps[j+*kk] - y) * (traps[j+*kk] - y);
		    
		    if ((d21<=truncate2) && (d22<=truncate2)) {
			p1p2 = exp(-(d21+d22) / 2 / *sigma / *sigma);
			if (wt[i] > 1)
			    p1p2 = 1 - pow(1-p1p2, wt[i]);   /* tentative 2012-12-15 */
		    }
		    else
			p1p2 = 0;
		    sump  += p1p2;
		    sumdp += p1p2 * sqrt(dij);
		}
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
  int    *wt,      /* integer trap weights */
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
	    if (wt[k] > 0) {
		d2k = (traps[k] - x) * (traps[k] - x) + (traps[k + *kk]-y) * (traps[k + *kk]-y);
		pk = pfn(*fn, d2k, 1.0, *sigma, *z, miscparm, truncate2);
		if (wt[k] > 1)
		    pk = 1 - pow(1-pk, wt[k]);
		sumd2k += pk * d2k;
		sumpk  += pk;
		pdot   *= (1 - pk);
	    }
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
  int    *detect,  /* scalar code 0 = multicatch, 1 = proximity */
  double *g0,      /* Parameter : detection magnitude */
  double *sigma,   /* Parameter : detection scale */
  int    *ss,      /* number of occasions */
  int    *kk,      /* number of traps */
  int    *mm,
  int    *wt,      /* integer trap weights */
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
    int varying = 0;

    if (*fn != 0)
    error ("invalid detection function in naivecap2");

    /* check to see if old code sufficient 2012-12-15 */
    for (k =0; k<*kk; k++) {
	if (wt[k] != *ss) {
	    varying = 1;
	    break;
	}
    }

    if (!varying) {
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
	    *value = 0;    
	else
	    *value = *ss * nsum / psum;
    }
    else {

/* abandon multicatch for now as requires more complex allowance for trap competition */
/* for multicatch need to accumulate Pr(caught) within each occasion, and 
   need entire usage for this */

	for (m=0; m<*mm; m++) {
	    pk = 0.0;
	    product = 1.0;
	    for (k=0; k<*kk; k++)
	    {
		if (wt[k] > 0) {
		    d2val = d2(m, k, mask, traps, *mm, *kk);
		    pk = *g0 * expmin(-d2val / 2 / *sigma / *sigma);
		    nsum += wt[k] * pk;  /* wt[k] opportunities, each probability pk */
		    if (wt[k] > 1)
			product *= pow(1 - pk, wt[k]); 
		    else
			product *= (1 - pk); 
		}
	    }
	    psum += 1 - product;   /* Pr capt animal at m over all traps & times */
	}
	if (psum<=0)
	    *value = 0;    /* failed */
	else
	    *value = nsum / psum;
    }
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
    size_t memrequest;  // "size_t is defined in stddef.h which the header defining R_alloc..."
    *resultcode = 1;

    memrequest = *nrow;
    memrequest *= *ncol;
    ytemp = (double *) R_alloc(memrequest, sizeof (double));

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

/* 
Code for pdf of one individual range centre, evaluated at an arbitrary set of points

The individual is selected with the first argument 'which'; this must be in range 1:n,
where n is the number of individuals used to fit the model.

Allows non-uniform D
Previously `pwuniform'

2014-08-04
*/

void fxIHP (
    int    *which,       /* which one: 1 <= which <= *nc */
    int    *xx,          /* number of new points at which f(X_i) requested */
    double *X,           /* new points at which f(X_i) requested */
    double *piX,         /* pimask at each requested point X */
    int    *like,        /* likelihood 0 full, 1 conditional */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
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
    double *dist2,       /* distances (optional: -1 if unused) */
    double *distX2,      /* distances to X points (optional: -1 if unused) */
    double *Tsk,         /* nk x s usage matrix */
    int    *markocc,     /* which are marking occasions? */
    int    *Tu,          /* detector x occasion matrices of no. sightings Tu */
    int    *Tm,          /* detector x occasion matrices of no. sightings Tm */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *pimask,      /* individual probability density; used for normalisation */
    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *PIA,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    double *miscparm,    /* miscellaneous parameters (cutval, normalization, etc.) */
    int    *normal,      /* code 0 don't normalise, 1 normalise */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard rate, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *minprob,     /* minimum value of P(detection history) */
    double *value,       /* return values */
    int    *resultcode   /* 0 if OK */
)
{
    int    i,m,x;
    int    *ng;       /* number per group */
    double *pmixg = NULL;
    double *pmixn = NULL;
    double temp;
    double sumprwi = 1.0;
    double prwi;

    double *gk = NULL;
    double *hk = NULL;
    double *gkx = NULL;
    double *hkx = NULL;
    double *detspec;
    double *detspecx;

//    double *hr = NULL;

    int    *hc0 = NULL;        
    int    *hindex = NULL;     
    double *h = NULL;          

    int    *hc0x;
    int    *hindexx;
    double *hx;            /* 2011-11-15 */ 
    int    *start = NULL;

    zfnptr zfn;
    prwfnptr prwfn;
    int cumk[maxnpoly];
    int nk = 0;
    int nd = 0;
    int nc1 = 0;
    int nv; 
    int ncol;
    int gpar = 2;    /* number of 'g' (detection) parameters */

    double scale = 1e6;

    double *pID;

    /*===============================================================*/

    /* MAINLINE */

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
    /*---------------------------------------------------------*/

    if (*gg > 1)
	error("fxi does not allow for groups");
    ng = (int *) R_alloc(*gg, sizeof(int));
    fillng(*nc, *gg, grp, ng);                                    /* number per group */

    nk = fillcumk(detect, *ss, kk, cumk);                          /* detections per polygon */

    /*---------------------------------------------------------*/

    /* Select detection function and detector */
    /* see utils.c for these functions */
    zfn = getzfn(*fn);                                          
    prwfn = getprwfn(detect, *ss);
    gpar = par3(*fn) + 2;

    /*---------------------------------------------------------*/

    /* Mark-resight */
    /* Update pID, gpar */
    pID = (double *)  S_alloc(*ss * *nmix, sizeof (double));
    if (*like == 1) ncol = nc1;
    else ncol = *gg;
    gpar = markresightini (*ss, *nmix, markocc, nk, ncol, PIA, *cc, gsbval,
		   pID, gpar);;

//    allsighting = (*like >= 5);
//    if ((Tu[0]>=0) || (Tm[0]>=0)) {
        /* storage for Pr(marked on or before s) */
//	pdots = (double *)  S_alloc(*ss * *nmix * *mm, sizeof (double));
//    }

    /*---------------------------------------------------------*/

    /* pimask has dual use: telemetry and allsighting */
    /*
    usepimask = (pimask[0] > -tol) && (!allsighting);
    if (usepimask || *like == 6) {
	pdot = (double *)  S_alloc(nc1 * *nmix * *mm, sizeof (double));
    }
    */
    /*---------------------------------------------------------*/

    /* 'start' is index to position of ancillary data for each
       individual, used for polygon, transect and sound detectors */
    if ( (detect[0] == 3) || (detect[0] == 4) ) 
	start = (int *) R_alloc(nc1 * *ss, sizeof(int));
    else if (((detect[0]>=5) && (detect[0]<=7)) || (detect[0]==12)) 
	start = (int *) R_alloc(nc1 * *ss * nk, sizeof(int));
    nd = getstart(detect, start, nc1, *nc, *ss, nk, w);
    /*---------------------------------------------------------*/

    pmixn = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    for (i=0; i < *gg * *nmix; i++) pmixg[i] = 1; /* default */
    for (i=0; i < nc1 * *nmix; i++) pmixn[i] = 1; /* default */
    getpmix(gpar, nc1, *nmix, knownclass, *nc, *cc, *ss, nk, grp, PIA, 
		   gsbval, pmixg, pmixn);
    nv = nval(detect[0], nc1, *cc, *ss, nk);

    /*---------------------------------------------------------*/
    /* distance options Sept 2014 */
    /* distances to mask points */
    /* skip for polygon or transect detectors */
    if (dist2[0] < 0) {
        if (anypolygon(detect,*ss) || anytransect(detect,*ss)) {
	    dist2 = (double *) S_alloc(1, sizeof(double));
	}
	else {
	    dist2 = (double *) S_alloc(*kk * *mm, sizeof(double));
	    makedist2 (*kk, *mm, traps, mask, dist2);
	}
    }
    else {
	squaredist(*kk, *mm, dist2);
    }

    /* distances to novel points X */
    /* skip for polygon or transect detectors */
    if (distX2[0] < 0) {
        if (anypolygon(detect,*ss) || anytransect(detect,*ss)) {
	    distX2 = (double *) S_alloc(1, sizeof(double));
	}
	else {
	    distX2 = (double *) S_alloc(*kk * *xx, sizeof(double));
	    makedist2 (*kk, *xx, traps, X, distX2);
	}
    }
    else {
	squaredist(*kk, *xx, distX2);
    }
    /*---------------------------------------------------------*/

    /* do this once */
    if (anyexclusive(detect,*ss)) {
        zerok (detect, *nc, *ss, nk, w);
    }

    /* optionally compute sumprwi, denominator used for normalisation */
    /* as at 2014-09-10 this does not allow hcov */
    if (*normal > 0) {
	gk = (double *) S_alloc(*cc * nk * *mm, sizeof(double)); /* S_alloc sets to zero */
	hk = (double *) S_alloc(*cc * nk * *mm, sizeof(double)); /* S_alloc sets to zero */
	detspec = (double *) R_alloc(nv, sizeof(double));
	precompute(detect, *fn, binomN, *ss, *kk, *mm, *cc, nk, cumk,  
		   traps, dist2, mask, gsbval, miscparm, detspec, gk, hk, 0); 
	/* space allocation and precomputation for exclusive detector types (traps) */
	if (anyexclusive(detect,*ss)) {
	    hc0 = (int *) R_alloc (*cc, sizeof(int));
	    hindex = (int *) S_alloc (nc1 * *ss, sizeof(int));
	    h = (double *) S_alloc (nc1 * *ss * *mm * *nmix, sizeof(double)); 
	    geth (nc1, *cc, *nmix, nk, *ss, *mm, PIA, hc0, hk, Tsk, h, hindex);
	} 
	getdetspec (detect, *fn, *nc, nc1, *cc, *nmix, nd, nk, *ss,  /* complete fill of detspec */
		    *kk, *mm, PIA, miscparm, start, detspec);
	
        sumprwi = 0;
        for (x = 0; x < *nmix; x++) {
            temp = 0;
            for (m = 0; m < *mm; m++) {
		prwi = prwfn (m, *which-1, x, w, xy, signal, PIA, gk, hk, binomN,
			      detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, zfn,
			      gsbval, traps, dist2, Tsk, mask, *minprob, pID);
		temp += prwi * pimask[m];
            }
            sumprwi += pmixg[x] * temp * scale;
        }    /* end of loop over mixtures */
        if (sumprwi < fuzz)
            error("zero prwi in external function fxIHP");
    }

    /*---------------------------------------------------------*/
    /* dynamically allocate memory for requested X locations   */
    /* use R_alloc for robust exit on interrupt                */
    /* S_alloc zeros memoy as well 2011-11-15                  */
    /* h for total hazard added 2011-11-15                     */

    gkx = (double *) S_alloc(*cc * nk * *xx, sizeof(double));
    hkx = (double *) S_alloc(*cc * nk * *xx, sizeof(double));
    detspecx = (double *) R_alloc(nv, sizeof(double));

    precompute(detect, *fn, binomN, *ss, *kk, *xx, *cc, nk, cumk,  
	       traps, distX2, X, gsbval, miscparm, detspecx, gkx, hkx, 0); 

    /* space allocation and precomputation for exclusive detector types (traps) */
    if (anyexclusive(detect,*ss)) {
        hc0x = (int *) R_alloc (*cc, sizeof(int));
        hindexx = (int *) S_alloc (nc1 * *ss, sizeof(int));
        hx = (double *) S_alloc (nc1 * *ss * *xx * *nmix, sizeof(double));
	geth (nc1, *cc, *nmix, nk, *ss, *xx, PIA, hc0x, hkx, Tsk, hx, hindexx);
    }
    else {
        hc0x = (int *) R_alloc (1, sizeof(int));
        hindexx = (int *) S_alloc (1, sizeof(int));
        hx = (double *) S_alloc (1, sizeof(double));
    }

    getdetspec (detect, *fn, *nc, nc1, *cc, *nmix, nd, nk, *ss,     /* complete filling of detspecx */
		*kk, *xx, PIA, miscparm, start, detspecx);


    R_CheckUserInterrupt();

    /* i indexes points at which to evaluate pdf */
    for (i=0; i< *xx; i++) {
        temp = 0;
	for (x=0; x<*nmix; x++) {
	    temp += 
		pmixg[x] * 
		prwfn (i, *which-1, x, w, xy, signal, PIA,
		       gkx, hkx, binomN, detspecx, hx, hindexx, *cc, *nc, nk, 
		       *ss, *xx, *nmix, zfn, gsbval, traps, distX2, Tsk, X, *minprob, pID) *
		piX[i] * scale;
	}
        value[i] = temp / sumprwi;
    }
    
    *resultcode = 0;   /* successful termination fxIHP */
}
/*==============================================================================*/

void chat(
    int    *like,        /* likelihood 0 full, 1 conditional etc. */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
    int    *distrib,     /* distribution of nc 0 Poisson, 1 binomial */
    int    *w,           /* capture histories (1:nc, 1:s, 1:k) */
    int    *grp,         /* group number for 0<=n<*nc   [full likelihood only] */
    int    *nc,          /* number of capture histories */
    int    *ss,          /* number of occasions */
    int    *kk,          /* number of traps */
    int    *mm,          /* number of points on mask */
    int    *gg,          /* number of groups */
    int    *nmix,        /* number of mixtures */
    int    *knownclass,  /* known membership of 'latent' classes; 1='unknown' */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *dist2,       /* distances (optional: -1 if unused) */
    double *Tsk,         /* nk x s usage matrix */
    int    *markocc,     /* which are marking occasions? */
    double *pimask,      /* pdf of marked animals like == 5,6 */
    double *mask,        /* x,y points on mask (first x, then y) */
    double *Dmask,       /* density at each point on mask, possibly x group */
    double *gsb0val,     /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    int    *cc0,         /* number of g0/sigma/b combinations for naive animals */
    int    *PIA0,        /* lookup which g0/sigma/b combination to use for given g, S, K */
    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameters (cutval, normalization, NT etc.) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    int    *nsim,        /* number of replicate simulations for chat */
    double *chat,        /* return vector chat(Tu), chat(Tm) */
    int    *resultcode   /* 0 if OK */
)

/* experimental code for standalone estimation of sighting c-hat by simulation */
/* mostly a cut-down version of secrloglik; not all inputs used */
{
    /* indices */
    int    i,m,x;

    /* numbers of individuals and detectors */
    int    nc1;
    int    cumk[maxnpoly];
    int    nk = 0;

    /* number of 'g' (detection) parameters */
    int    gpar;    

    /* number of marking occasions and related */
    double *pID;
    double *pdots = NULL;
    int    allsighting; 
    double a0[maxnmix];

    /* pre-computed detection probabilities */
    double *gk0 = NULL;
    double *hk0 = NULL;

    /* mixture membership probability */
    double *pmixg = NULL;
    double *pmixn = NULL;

    double *detspec = NULL;
    int    nv;
    int    ncol;
    /*    double pdt; */

    /*-------------------------------------------------------*/

    /* MAINLINE */

    *resultcode = 1;

    /*---------------------------------------------------------*/
    /* Adjust for zero animals                                 */
    /* (otherwise fails with nc = 0 because *nc also used to   */
    /* dimension arrays etc. so replace nc1 = max(1, *nc)      */
    if (*nc) nc1 = *nc; else nc1 = 1;

    /*---------------------------------------------------------*/

    /* Detections per polygon */
    nk = fillcumk(detect, *ss, kk, cumk);                            
    /*---------------------------------------------------------*/

    gpar = par3(*fn) + 2;

    /* Mark-resight */
    /* Update pID, gpar, storage for Pr(marked on or before s) */
    pID = (double *)  R_alloc(*ss * *nmix, sizeof (double));
    pdots = (double *)  S_alloc(*ss * *nmix * *mm, sizeof (double));
    if (*like == 1) ncol = nc1;
    else ncol = *gg;
    gpar = markresightini (*ss, *nmix, markocc, nk, ncol, PIA0, *cc0, 
			   gsb0val, pID, gpar);
    allsighting = (*like >= 5);

    /*---------------------------------------------------------*/

    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    pmixn = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    for (i=0; i < (*gg * *nmix); i++) pmixg[i] = 1; /* default */
    for (i=0; i < (nc1 * *nmix); i++) pmixn[i] = 1; /* default */
    getpmix(gpar, nc1, *nmix, knownclass, *nc, *cc0, *ss, nk, 
		   grp, PIA0, gsb0val, pmixg, pmixn);

    nv = nval(detect[0], nc1, *cc0, *ss, nk);
    detspec = (double *) R_alloc(nv, sizeof(double));
    gk0 = (double *) S_alloc(*cc0 * nk * *mm, sizeof(double));
    hk0 = (double *) S_alloc(*cc0 * nk * *mm, sizeof(double));

    /*-----------------------------------------------------------*/
    if (dist2[0] < 0) {
        if (anypolygon(detect, *ss) || anytransect(detect,*ss)) {
	    dist2 = (double *) S_alloc(1, sizeof(double));
	}
	else {
	    dist2 = (double *) S_alloc(*kk * *mm, sizeof(double));
	    makedist2 (*kk, *mm, traps, mask, dist2);
	}
    }
    else {
	squaredist(*kk, *mm, dist2);
    }
    /*-----------------------------------------------------------*/

    precompute(detect, *fn, binomN, *ss, *kk, *mm, *cc0, nk, cumk,   
	       traps, dist2, mask, gsb0val, miscparm, detspec, gk0, hk0, 0); 

    /*
     if (anyexclusive(detect,*ss)) 
        error ("exclusive detectors not supported");
    */
    
    if ((detect[0]==5) || (detect[0]>7))
        error ("detector type not supported for chat");
    if (*nsim < 2)
	error ("nsim for chat must be at least 2, and preferably much more!");

    // no groups
    *resultcode = 0;
    if (*like < 5) { /* otherwise getpdots not used by sightingchat */
	for (x=0; x < *nmix; x++) {
	    for (m=0; m < *mm; m++)  {
		getpdots (m, 0, markocc, x, ncol, PIA0, gk0, hk0, detect,
			binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, pdots);
	    }
	}
    }
    else
	*resultcode = filla0 (*like, 0, markocc, ncol, PIA0, gk0, hk0, detect, binomN, Tsk, 
			      *ss, nk, *mm, *cc0, *nmix, gsb0val, allsighting, pimask, 
			      *area, a0);
    if (*resultcode>0) {
	// Rprintf("bad a0\n");
	return;
    }

    *resultcode = sightingchat (*like, detect, binomN, *nc, *ss, nk, *cc0, *nmix, pmixg, 
				*mm, Dmask, pimask, *area, pID, markocc, *nc, ncol, PIA0,
				gk0, hk0, Tsk, a0, *distrib, *nsim, pdots, chat);
}
/*==============================================================================*/
