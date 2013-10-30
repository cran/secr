/*
   External procedures for secr package

   can compile with gcc 4.6.3 :
   gcc -Ic:/R/R-2.15.2/include -c secr.c -Wall -pedantic -std=gnu99

   [confirmed 2012-11-13 0640]

*/
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
/* 2013-07-02 slight adjustments for param=2 (treated as param=0) */

/*
        *detect may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  exclusive polygon detector
        4  exclusive transect detector
        5  signal detector
        6  polygon detector
        7  transect detector
        8  times  (undocumented)
        9  cue    (undocumented)
	12 signalnoise
*/

#include "secr.h"
#include <time.h>

/*==============================================================================*/

void R_CheckUserInterrupt(void);

/*==============================================================================*/

double pndot (int m, int n, int s1, int s2, int x, int ncol, int PIA0[],
	      double gk0[], int detect, int binomN, double Tsk[], int ss, int kk, 
	      int cc0, int nmix, double gsb0val[], int param)

/*
    probability animal at point m on mask is caught
    n may indicate group (full likelihood; ncol= number of groups) or
    individual (conditional likelihood; ncol= number of individuals)
    aligned with secrloglik 2009 06 25

    2009 10 24 adjusted to allow summation over qq < ss
    2009 11 12 'kk' should be number of parts for polygon detectors
    2011 01 04 'param = 1' for GR parameterisation of multi
    2012 12 17 Tsk
    2012 12 25 for x=0, dbinom(x, N, 1-(1-p)^Tsk) = dbinom(x, N * Tsk, p)
    2012 12 25 DOES NOT ALLOW binomN < 0
*/
{
    int k,s,c;
    int gi;
    int wxi;
    double pp;
    double p0;
    double Ei;
    double g1;
    double Tski;
    if (binomN < 0) error ("negative binomN not allowed in C fn pndot");
    pp = 1;
    if (param == 1) {
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
                    Ei += (1-gk0[gi]) / p0;
                }
            }
            pp *= (1 - p0 * exp(-1/Ei));
        }
        return (1 - pp);
    }
    else {
        for (s=s1-1; s<s2; s++) {
            for (k=0; k< kk; k++) {
                wxi = i4(n,s,k,x,ncol,ss,kk);
                c = PIA0[wxi] - 1;
                if (c >= 0) {    /* drops unset traps */
                    gi = i3(c,k,m,cc0,kk);
		    g1 = gk0[gi];
		    Tski = Tsk[s * kk + k];
		    /* expect binomN = 1 if not count detector */
		    if (fabs(Tski-1) > 1e-10) {                 /* effort <> 1.0 */
			if ((detect < 8) & (detect != 5))  {
			    if (binomN == 0)
				g1 = exp(-Tski * g1);           /* Poisson */
			    else
				g1 = pow(1-g1, Tski * binomN);  /* Binomial */
			}
			else error("no effort adjustment for detector type");
		    }
		    else {
			if (binomN == 0)
			    g1 = exp(-g1);              /* Poisson */
			else  if (binomN == 1)
			    g1 = 1 - g1;                /* Bernoulli */
			else {
			    g1 = pow(1-g1, binomN);     /* Binomial */
			}
		    }
                    pp *= g1;
                }
            }
        }
        return (1 - pp);
    }
}
/*===============================================================*/

double prwimulti
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

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
    double Tski;
    for (s=s1-1; s<s2; s++)
    {
        htemp = h[i3(x,m,hindex[s*nc + n],nmix, mm)];
        if (htemp < fuzz) { result = 0; break; }
        k = w[nc * s + n];
        if (k < 0) {dead=1; k=-k;}  /*  1 <= k <= kk */
        if (k > 0) {
	    Tski = Tsk[s * kk + k-1];
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = -Tski * log (1 - gk[gi]);
            pks *= (1-expmin(-htemp)) / htemp;
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
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)
/*
    Probability of capture history n (0 <= n < *nc)
    given that animal's range centre is at m
    MULTI-CATCH DETECTOR, GARDNER & ROYLE PARAMETERISATION
    NO ALLOWANCE FOR EFFORT Tsk 2012-12-17
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

double prwiprox
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    BINARY PROXIMITY DETECTOR distinguished from count detector 2012-12-22
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of trap      0 <= k < kk  */
    int c, gi, wi, wxi;
    int count;
    int dead = 0;
    double result = 1.0;
    double g1;
    double Tski;

    for (s=s1-1; s<s2; s++) {
        for (k=0; k<kk; k++) {
            wi = i3(n, s, k, nc, ss);
            wxi = i4(n, s, k, x, nc, ss, kk);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            c = PIA[wxi] - 1;
            if (c >= 0) {                                       /* skip if this trap not set */
                gi  = i3(c, k, m, cc, kk);
		g1 = gk[gi];
		Tski = Tsk[s * kk + k];
		if (fabs(Tski-1) > 1e-10) /* not unity */
		    g1 = 1 - pow(1 - g1, Tski);
		if (count)                                      /* Bernoulli */
		    result *= g1;
		else 
		    result *= (1-g1);
                if (result < minp) {result = minp; break;}      /* truncate */
            }
        }
        if (result <= minp) {result = minp; break;}             /* truncate */
        if (dead==1) break;
    }
    return (result);
}
/*=============================================================*/

double prwicount
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

/*
    Probability of capture history n (0 <= n < nc)
    given that animal's range centre is at m
    COUNT DETECTOR, excludes BINARY PROXIMITY
*/
{
    int s;   /* index of occasion  0 <= s < ss  */
    int k;   /* index of trap      0 <= k < kk  */
    int c, gi, wi, wxi;
    int count;
    int dead = 0;
    double result = 1.0;
    double g1;
    double Tski;

    if (binomN < 0) error ("negative binomN not allowed in C fn prwicount");
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<kk; k++) {
            wi = i3(n, s, k, nc, ss);
            wxi = i4(n, s, k, x, nc, ss, kk);
            count = w[wi];
            if (count<0) {dead = 1; count = -count;}
            c = PIA[wxi] - 1;
            if (c >= 0) {                               /* skip if this trap not set */
                gi  = i3(c, k, m, cc, kk);
		g1 = gk[gi];
		Tski = Tsk[s * kk + k];
		if (binomN == 0) {                      /* Poisson */
		    result *= dpois(count, Tski * g1, 0); 
		}
		else if (binomN == 1) {                 /* Binomial, size from Tsk */ 
		    result *= countp (count, round(Tski), g1);
		}
		else if (binomN > 1) {                  /* Binomial, specified size */
		    if (fabs(Tski-1) > 1e-10) /* not unity */
			g1 = 1 - pow(1 - g1, Tski);
		    result *= countp (count, binomN, g1);
		}
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
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

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
		    g1 = gk[gi];
                    result *= countp (0, binomN, g1);
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
			    g1 = gk[gi];
			    result *= countp (1, binomN, g1);
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
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

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
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

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
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)
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
                lambda = Tsk[s * kk + k] * gk[gi];

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
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

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
    double g1;
    double Tski;

    if (binomN < 0) error ("negative binomN < 0 not allowed in C fn prwipoly");

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
		g1 = gk[gi];
		Tski = Tsk[s * kk + k];
		if (binomN == 0) {                      /* Poisson */
		    result *= dpois (count, Tski * g1, 0);
		}
		else if (binomN == 1) {                 /* Binomial, size from Tsk */ 
		    result *= countp (count, round(Tski), g1);
		}
		else {                                  /* Binomial */ 
		    if (fabs(Tski-1) > 1e-10) 
			g1 = 1 - pow(1 - g1, Tski);
		    result *= countp (count, binomN, g1);
		}
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
    /* Rprintf("prwipolygon %12.9f\n", result); */
    return (result);
}
/*=============================================================*/

double prwipolygonX
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)

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
    double Tski;

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
	    Tski = Tsk[s * kk + k-1];
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = -Tski * log (1 - gk[gi]);
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

double prwitransect
    (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
     int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
     int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
     double traps[], double Tsk[], double mask[], double minp)

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
    double g1;
    double Tski;

    nk = round(detspec[0]);
    nd = round(detspec[1]);
    for (s=s1-1; s<s2; s++) {
        for (k=0; k<nk; k++) {
            wi = i3(n,s,k,nc,ss);
            count = w[wi];
            if (count<0) {dead=1; count=-count;}
            wxi = i4(n,s,k,x,nc,ss,kk);
            c = PIA[wxi] - 1;
            if (c >= 0) {                               /* skip if this transect not used */
		Tski = Tsk[s * kk + k];
                gi  = i3(c,k,m,cc,nk);
                g1 = gk[gi];
		if (binomN == 0) {                      /* Poisson */
		    result *= countp (count, 0, Tski * g1);
		}
		else if (binomN == 1) {                 /* Binomial, size from Tsk */ 
		    result *= countp (count, round(Tski), g1);
		}
		else {                                  /* Binomial, size from binomN */ 
		    if (fabs(Tski-1) > 1e-10) 
			g1 = 1 - pow(1 - g1, Tski);
		    result *= countp (count, binomN, g1);
		}
                if ((result > minp) && (count>0)) {     /* avoid underflow 2009 11 16 */
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

double prwitransectX
   (int m, int n, int s1, int s2, int x, int w[], double xy[], double signal[],
    int PIA[], double gk[], int binomN, double detspec[], double h[], int hindex[],
    int cc, int nc, int kk, int ss, int mm, int nmix, gfnptr gfn, double gsbval[],
    double traps[], double Tsk[], double mask[], double minp)
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
    double Tski;

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
	    Tski = Tsk[s * kk + k-1];
            c = PIA[i4(n, s, k-1, x, nc, ss, kk)] - 1;
            gi = i3(c, k-1, m, cc, kk);
            pks = -Tski * log (1 - gk[gi]);
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

void pdotpoint (double *xy, int *nxy, double *traps, int *detect, 
                double *Tsk, int *kk, int *fn, double *par, int *nocc, 
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

    if (*fn>18)
        if (*fn != 20) error("pdotpoint requires detectfn < 18");
    g0 = par[0];
    sigma = par[1];
    if (!((*fn == 0) || (*fn == 2) || (*fn == 4) || (*fn == 9) || (*fn == 20)
	  || (*fn == 14) || (*fn == 16)))
        z = par[2];
    if ((*fn == 10) || (*fn == 11) || (*fn == 12) || (*fn == 13))
        cutval[0] = par[3];

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
    else for (i=0; i<*nxy; i++) {
	tempval = 1;
	for (s=0; s<*nocc; s++) {
	    for (k=0; k<*kk; k++) {
		Tski = Tsk[s * *kk + k];
		if (Tski > 1e-10) {
		    dk2 = (xy[i]-traps[k]) * (xy[i]-traps[k]) +
			(xy[i + *nxy]-traps[k+ *kk]) * (xy[i + *nxy]-traps[k+ *kk]);
		    p = pfn(*fn, dk2, g0, sigma, z, cutval, *w2);

		    if (*detect == 2) {    /* counts */
			if (*binomN == 0)
			    p = 1 - countp(0, 0, Tski * p);
			else if (*binomN == 1)
			    p = 1 - countp(0, round(Tski), p);
			else {
			    if (fabs(Tski-1) > 1e-10)
				p = 1 - pow(1-p, Tski);
			    p = 1 - countp(0, *binomN, p);
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
	value[i] = 1 - tempval;
    }
}

/*=============================================================*/

void pdotpoly (double *xy, int *nxy, double *traps, int *detect, double *Tsk, int *nk,
    int *kk, int *fn, double *par, int *nocc, int *binomN, double *value)
/* ignores effort except 0/1 */
/* incorporated pdottransect 2012-12-24 */
/* NEED TO CLARIFY HAZARD FORMULATION 2012-12-24 */
{
    int i,k,s;
    double hk = 0;
    double sumhk;
    double tempval;
    int *cumk;
    double stdint = 1.0;
    double *ex;
    double Tski = 1.0;
    ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));

    cumk = (int *) R_alloc(*nk+1, sizeof(int));
    cumk[0] = 0;
    for (i=0; i<*nk; i++)
        cumk[i+1] = cumk[i] + kk[i];

    if ((*detect == 3) | (*detect == 6))          /* polygon */
	stdint = gintegral(*fn, par);
    else if ((*detect == 4) | (*detect == 7))     /* transect */
	stdint = gintegral1(*fn, par);
    else 
	error ("unrecognised detector type in pdotpoly");

    for (i=0; i<*nxy; i++) {
	tempval = 1.0;
	for (s=0; s<*nocc; s++) {
	    sumhk = 0.0;
	    for (k=0; k<*nk; k++) {
                Tski = Tsk[s * *nk + k];
		if (Tski > 1e-10) {
		    if ((*detect == 3) | (*detect == 6))
			hk = par[0] * integral2D (*fn, i, 0, par, 1, traps, xy,
			    cumk[k], cumk[k+1]-1, cumk[*nk], *nxy, ex) / stdint;
		    else if ((*detect == 4) | (*detect == 7))
			hk = par[0] * integral1D (*fn, i, 0, par, 1, traps, xy,
			    cumk[k], cumk[k+1]-1, cumk[*nk], *nxy, ex) / stdint;
		    sumhk += hk;
		}
	    }
	    tempval *= countp(0, *binomN, sumhk);
	}
	value[i] = 1 - tempval;
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

void precompute(
    int *detect, 
    int *fn, 
    int *binomN,
    int *kk,
    int *mm,
    int *cc,
    int nk,
    int cumk[],
    double traps[],
    double mask[],
    double gsbval[],
    double miscparm[],
    double detspec[],
    double gk[]
 ) {
    /*---------------------------------------------------------*/
    /* populate pre-computed gk and gk0 arrays                            */
    /*
        *detect may take values -
        0  multi-catch traps
        1  binary proximity detectors
        2  count  proximity detectors
        3  exclusive polygon detector
        4  exclusive transect detector
        5  signal detector
        6  polygon detector
        7  transect detector
        8  times  (undocumented)
        9  cue    (undocumented)
	12 signalnoise
    */
    /*---------------------------------------------------------*/

    int c,m,k,gi;

    double par[4];   /* passing parameter values to integr fn  */
    double stdint = 1;
    gfnptr gfn;
    double *ex; 

    gfn = getgfn(*fn);

    /* all these detectors are Bernoulli, or similar           */
    /* so we override binomN                                   */
    if ((*detect==0) || (*detect==1) || (*detect==3) || 
        (*detect==4) || (*detect==5) || (*detect==12)) {
        *binomN = 1;
    }

    if (*detect == 0) {
        for (k=0; k<nk; k++) {
            for (m=0; m<*mm; m++) {
                for (c=0; c<*cc; c++) {
		    gi = i3(c,k,m,*cc,nk);
                    gk[gi] = gfn(k, m, c, gsbval, *cc, traps, 
				 mask, nk, *mm, miscparm);
                }
            }
        }
    }
    else if ((*detect == 1) || (*detect == 2) || (*detect==5) || (*detect==8)
	     || (*detect==9) || (*detect==12)) {
        for (k=0; k<nk; k++) {
            for (m=0; m<*mm; m++) {
                for (c=0; c<*cc; c++) {
                    gi = i3(c,k,m,*cc,nk);
		    gk[gi] = gfn(k, m, c, gsbval, *cc, traps, mask, nk, *mm, miscparm);
                }
            }
        }
    }
    else if ((*detect == 3) || (*detect == 6)) {
	ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
        for (c=0; c<*cc; c++) {
            par[0] = gsbval[c];
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];
            stdint = gintegral(*fn, par);
            detspec[2+c] = stdint;               /* passed to prwipolygon */
            for (k=0; k<nk; k++) {               /* over parts */
                for (m=0; m<*mm; m++) {
                    gi = i3(c, k, m, *cc, nk);
			gk[gi] = par[0] * integral2D (*fn, m, 0, par,
                            1, traps, mask, cumk[k], cumk[k+1]-1, cumk[nk],
						      *mm, ex) / stdint;
                }
            }
        }
    }
    else if ((*detect == 4) || (*detect == 7)) {
	ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
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
    		          traps, mask, cumk[k], cumk[k+1]-1, cumk[nk], *mm, ex) /
                            stdint;
                }
            }
        }
    }
    else error ("unrecognised detector type in external C fn secrloglik");
} /* end precompute */
/*=============================================================*/

prwfnptr getprwfn (int detect, int param) 
{
    prwfnptr prwfn;
    prwfn = prwicount; /* default */
    if ((detect == 0) && (param == 1))
        prwfn = prwimultiGR;
    else if (detect == 0)
        prwfn = prwimulti;
    else if (detect == 1)
        prwfn = prwiprox;
    else if (detect == 3)
        prwfn = prwipolygonX;
    else if (detect == 4)
        prwfn = prwitransectX;
    else if ((detect == 5) || (detect==9)) 
	prwfn = prwisignal;
    else if (detect == 6)
        prwfn = prwipolygon;
    else if (detect == 7)
        prwfn = prwitransect;
    else if (detect == 8)
        prwfn = prwitimes;
    else if (detect == 12) {
	prwfn = prwisignalnoise;   /* experimental 2012-02-07 */
    }
    return(prwfn);
} /* end getprwfn */
/*=============================================================*/

int getstart(int *detect, int start[], int nc1, int *nc, 
	     int *ss, int nk, int w[]) {

    /*---------------------------------------------------------*/
    /* Identify start positions of ancillary data for each     */
    /* animal                                                  */

    int nd = 0;
    int i,s,k,wi;

    /* max one detection per occasion */
    if ( (*detect == 3) || (*detect==4) ) {
        /* start[z] indexes the row in xy
           for each detection z, where z is w-order (is) */
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
    return(nd);
}
/*=============================================================*/

void getdetspec (int *detect, int *fn, int *nc,  int nc1, 
		 int *cc, int *nmix, int nd, int nk, int *ss, int *kk, int *mm, 
		 int PIA[], double miscparm[], int start[], double detspec[]) {

    /* detector-specific data passed later to prwi functions */

    /* mixtures are group-specific for full likelihood, and     */
    /* individual-specific for conditional likelihood           */

    int i;

    if ((*detect == 3) || (*detect == 4)) {
        detspec[0] = (double) nk;
        detspec[1] = (double) nd;
        /* pass index of first detection of each animal + occasion */
        /* maximum 1 per polygon as exclusive */
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
        detspec[1] = (double) nd;                     /* number of detections or detectors??*/
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
        /* pass index of first detection of each animal + occasion + polygon */
        /* possibly more than one per polygon */
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[2+*cc+i] = (double) start[i];
    }
    else if (*detect == 8) {
        for (i=0; i< (*nc * *ss * nk); i++)
            detspec[i] = (double) start[i];
    }
}
/*=============================================================*/

void geth (int nc1, int *cc, int *nmix, int nk, int *ss, int *mm, 
	   int PIA[], int hc0[], double gk[], double Tsk[],
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
        for (s=0; s < *ss; s++) {
	    PIAval0 = PIA[i4(n,s,0,0, nc1, *ss, nk)];
	    for (k=1; k<nk; k++) {
                PIAvalk = PIA[i4(n,s,k,0, nc1, *ss, nk)];
		if (PIAval0 < 0)
		    PIAval0 = PIAvalk;
		else if (PIAvalk>0) {
		    if (PIAval0 != PIAvalk) {
			fullns = 1;
			break;
		    } 
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
		for (k=0; k < nk; k++) {
		    Tski = Tsk[s * nk + k];
		    for (x = 0; x < *nmix; x++) {
			c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
			for (m = 0; m < *mm; m++) { 
			    if (c >= 0) {
				gi = i3(c,k,m,*cc,nk);
				h[i3(x,m,hi,*nmix, *mm)] += Tski * hazard(gk[gi]);
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
		    for (k=0; k < nk; k++) {
			Tski = Tsk[s * nk + k];
			for (x = 0; x < *nmix; x++) {
			    c = PIA[i4(n,s,k,x, nc1, *ss, nk)]-1; 
			    for (m=0; m< *mm; m++) { 
				if (c >= 0) {
				    gi = i3(c,k,m,*cc,nk);
				    h[i3(x,m, hc0[c0],*nmix, *mm)] += Tski * hazard(gk[gi]);
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

int getpmix(int gpar, int nc1, int *nmix, int knownclass[], int *nc, int *cc, int *ss, 
            int nk, int grp[], int PIA[], double gsbval[], double pmixg[], 
            double pmixn[]){

    /*-------------------------------*/
    /* mixture proportions           */
    /* by group and by animal        */

    int g, n, c, x, wxi;
    double pmix;
    if (*nmix>1) {
        /* one extra real parameter for h2 */
        gpar++;
        for (n=0; n<nc1; n++) {
            for (x=0; x<*nmix; x++) {
                wxi = i4(n,0,0,x,*nc,*ss,nk);
                c = PIA[wxi] - 1;
		pmix = gsbval[*cc * (gpar-1) + c];
		g = grp[n]-1;
		pmixg[*nmix * g + x] = pmix;
                
		if (knownclass[n] > 1) {
		    if (knownclass[n] == (x+2))   /* knownclass=2 maps to x=0 */
			pmixn[*nmix * n + x] = 1;
		    else 
			pmixn[*nmix * n + x] = 0;
		}
		else
		    pmixn[*nmix * n + x] = pmix;
            }
        }
    }
    return(gpar);  /* incremented by 1 if mixture model */
}
/*=============================================================*/

int nval(int *detect, int nc1, int *cc, int *ss, int nk) {
    /* compute and allocate the space needed for ancillary data */
    int nval;
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
    return(nval);
}
/*=============================================================*/

int fillcumk(int *detect, int kk[], int cumk[]){
    /* determine number of polygons if polygon detector */
    /* for polygon detectors, kk is vector ending in zero */
    /* and first nk+1 elements of vector cumk are filled */
    int i;
    int nk = 0;
    if ((*detect==3) || (*detect==4) || (*detect==6) || (*detect==7)) {
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
        nk = kk[0];

    return(nk);
}
/*=============================================================*/

void fillng(int nc, int gg, int *grp, int ng[]) {
    int g, n;
    /* Count number per group (not used for CL)                */
    /* Assume histories sorted by group = individual           */
    /* CH are numbered 0 <= n < *nc in C code                  */
    for (g=0; g<gg; g++)
        ng[g] = 0;
    for (n=0; n<nc; n++) { 
        g = grp[n] - 1; 
        ng[g]++;
    }
}
/*=============================================================*/

void integralprw1 (
    int    *detect,    /* detector 0 multi, 1 proximity etc. */
    int    *param,     /* parameterisation 0 Borchers & Efford 1 Gardner & Royle */
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
    double *Tsk,       /* nk x s usage matrix */
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
    int i,n,m,x;
    double asum = 0;
    double *gk0;
    double D = 1.0;
    int cumk[maxnpoly];
    int nk = 0;
    int nc1;
    double *pmixg;      /* proportion in each mixture by group*/
    double *pmixn;      /* proportion in each mixture by individual*/
    int gpar = 2;
    int    *ng;       /* number per group */
    double *detspec = NULL;
    int nv;

    *resultcode = 1;                   /* generic failure */

    /* groups are not used in functions that call integralprw1 */
    /* but we go through the motions...*/
    ng = (int *) R_alloc(*gg, sizeof(int));
    fillng(*nc, *gg, grp, ng);                                    /* number per group */
    nk = fillcumk(detect, kk, cumk);                            /* detections per polygon */

    /*---------------------------------------------------------*/
    /* Adjust for zero animals                                 */
    /* (otherwise fails with nc = 0 because *nc also used to   */
    /* dimension arrays etc. so replace nc1 = max(1, *nc)      */
    if (*nc == 0) 
        nc1 = 1;
    else
        nc1 = *nc;
    /*---------------------------------------------------------*/

    nv = nval(detect, nc1, cc0, ss, nk);

    detspec = (double *) R_alloc(nv, sizeof(double));

    /* Allocate space for array of naive detection probability */
    gk0 = (double *) R_alloc(*cc0 * nk * *mm, sizeof (double));

    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    pmixn = (double *)  R_alloc(*nc * *nmix, sizeof (double));
    for (i=0; i < *gg * *nmix; i++) pmixg[i] = 1; /* default */
    for (i=0; i < nc1 * *nmix; i++) pmixn[i] = 1; /* default */
    gpar = getpmix(gpar, nc1, nmix, knownclass, nc, cc0, ss, nk, 
		   grp, PIA0, gsb0val, pmixg, pmixn);
    precompute(detect, fn, binomN, kk, mm, cc0, nk, cumk,   
	       traps, mask, gsb0val, miscparm, detspec, gk0); 

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
                    asum += pmixn[*nmix * n + x] * D * pndot (m, n, 1, *ss, x,
   		        *ncol, PIA0, gk0, *detect, *binomN, Tsk, *ss, nk, *cc0, 
			*nmix, gsb0val, *param);
                }
            }
        }
        a[n] = *area * asum;
    }
    *resultcode = 0;                   /* successful completion */
}

/*==============================================================================*/

clock_t timestamp(clock_t ticks1, int *counter) {
    clock_t ticks2;
    ticks2 = clock();
    Rprintf("check %2d: time used %12d ticks\n", *counter, ticks2-ticks1); 
    *counter = *counter + 1;
    return(ticks2);
}

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
    int    *knownclass,  /* known membership of 'latent' classes; 1='unknown' */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *Tsk,         /* nk x s usage matrix */
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
    /* indices */
    int    i,n,g,k,m,s,x;
    int    wxi;

    /* miscellaneous */
    double temp, tempg, tempsum, tempp, templog, prwi;

    /* group arrays */
    int    *ng;       /* number per group */
    int    *nm;       /* number per known mixture class */
    double *sumD;     
    double *sumDp;    

    /* generalised detection and detector functions */
    gfnptr gfn;
    prwfnptr prwfn;

    /* numbers of individuals and detectors */
    int    nc1;
    int    cumk[maxnpoly];
    int    nk = 0;

    /* number of 'g' (detection) parameters */
    int    gpar = 2;    

    /* number of detections */
    int    nd = 0;  

    /* pre-computed detection probabilities */
    double *gk = NULL;
    double *gk0 = NULL;

    /* mixture membership probability */
    double *pmixg = NULL;
    double *pmixn = NULL;
    double pmixnx;
    int known = 0;

    /* passing data to prwfn */
    double *detspec = NULL;
    int    *start = NULL;

    /* total hazard computation 2011-11-15*/
    int    *hc0 = NULL;        
    int    *hindex = NULL;     
    double *h = NULL;          

    /* CL only */
    int    indiv = 0;     /* indicate detection varying between individuals */
    double asum[maxnmix];
    double pd;
    int    nested = 0;  /* 0 = not nested, 1 = nested (CL only) */
    int    cumng = 0;   /* used for 'nested' */
    double dp;          /* Poisson density (cuerate) */


    double pdt = 0;
    int    nv;
    double tempN;

    /*-------------------------------------------------------*/

    /* MAINLINE */

     clock_t ticks;
     int timing = 0;
     int counter = 0;
     ticks = clock();

     /* 0 */
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

    /* Under development 2011-01-11                            */
    /* *gg > 1 with CL is used for clustered model             */
    if ((*gg > 1) && (*like == 1)) nested = 1;                 
    /*---------------------------------------------------------*/

    /* Number per known mixture, and total in known class */
    /* On input, knownclass 1 is used for 'unknown' : */
    /* knownclass 1 -> nm[0] 'unknown' */
    /* knownclass 2 -> nm[1] 'latent class 1' */
    /* knownclass 3 -> nm[2] 'latent class 2' */
    nm = (int *) R_alloc(*nmix + 1, sizeof(int));
    fillng(*nc, *nmix+1, knownclass, nm);                           
    if (*nmix > 1) {
	for (x = 1; x < (*nmix+1); x++)
	known += nm[x];
    }

    /* Number per group */
    ng = (int *) R_alloc(*gg, sizeof(int));
    fillng(*nc, *gg, grp, ng);                                    

    /* Detections per polygon */
    nk = fillcumk(detect, kk, cumk);                            

    /* Select detection function and detector */
    /* see utils.c for these functions */
    gfn = getgfn(*fn);                                          
    prwfn = getprwfn(*detect, *param);

    /* 'start' is index to position of ancillary data for each
       individual, used for polygon, transect and sound detectors */
    if ( (*detect == 3) || (*detect==4) ) 
	start = (int *) R_alloc(nc1 * *ss, sizeof(int));
    else if (((*detect>=5) && (*detect<=9)) || (*detect==12)) 
	start = (int *) R_alloc(nc1 * *ss * nk, sizeof(int));
    nd = getstart(detect, start, nc1, nc, ss, nk, w);

    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    pmixn = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    for (i=0; i < (*gg * *nmix); i++) pmixg[i] = 1; /* default */
    for (i=0; i < (nc1 * *nmix); i++) pmixn[i] = 1; /* default */
    gpar = getpmix(gpar, nc1, nmix, knownclass, nc, cc, ss, nk, 
		   grp, PIA, gsbval, pmixg, pmixn);

/* DEBUG 
    for (g=0; g < *gg; g++)
    for (x=0; x < *nmix; x++)
    Rprintf("g %d x %d pmixg %6.4f \n", g, x, pmixg[*nmix * g + x]);
    for (n=0; n < *nc; n++)
    for (x=0; x < *nmix; x++)
    Rprintf("n %d x %d pmixn %6.4f \n", n, x,  pmixn[*nmix * n + x]);
end DEBUG */

    nv = nval(detect, nc1, cc, ss, nk);
    detspec = (double *) R_alloc(nv, sizeof(double));

    gk = (double *) S_alloc(*cc * nk * *mm, sizeof(double));    /* S_alloc sets to zero */
    gk0 = (double *) S_alloc(*cc0 * nk * *mm, sizeof(double));

    /* 1 */
    if (timing) ticks = timestamp(ticks, &counter);

    precompute(detect, fn, binomN, kk, mm, cc, nk, cumk,    
	       traps, mask, gsbval, miscparm, detspec, gk); 
    precompute(detect, fn, binomN, kk, mm, cc0, nk, cumk,   
	       traps, mask, gsb0val, miscparm, detspec, gk0); 

    /* 2 */
     if (timing) ticks = timestamp(ticks, &counter);

    R_CheckUserInterrupt();
    if ((*detect==0) || (*detect==3) || (*detect==4)) {   /* exclusive detectors require hazard */ 
        hc0 = (int *) R_alloc (*cc, sizeof(int));
        hindex = (int *) S_alloc (nc1 * *ss, sizeof(int));
	h = (double *) S_alloc (nc1 * *ss * *mm * *nmix, sizeof(double)); 
	geth (nc1, cc, nmix, nk, ss, mm, PIA, hc0, gk, Tsk, h, hindex);
    } 

    getdetspec (detect, fn, nc, nc1, cc, nmix, nd, nk, ss,     /* complete filling of detspec */
		kk, mm, PIA, miscparm, start, detspec);

    /*-------------------------*/
    /* Now evaluate likelihood */
    /*-------------------------*/

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
	    /* Rprintf("all same\n"); */
            /* all individuals the same */
            /* save time by doing this once, rather than inside n loop */
            for (x=0; x<*nmix; x++) {
                asum[x] = 0;
                for (m=0; m<*mm; m++) {
                    if (nested) {
                        pd = pndot (m, 0, 1, *ss, x, *nc, PIA0, gk0, *detect, *binomN, Tsk,
                            *ss, nk, *cc0, *nmix, gsb0val, *param);
                        asum[x] += pndotgrp (pd, miscparm[0], 0.0001);
                    }
                    else {
			asum[x] += pndot (m, 0, 1, *ss, x, *nc, PIA0, gk0, *detect, *binomN, 
                           Tsk, *ss, nk, *cc0, *nmix, gsb0val, *param);
                    }
                }
            }
        }
        /* else asum calculated for each individual in loop below */
        *value = 0;

    /* 3 */
     if (timing) ticks = timestamp(ticks, &counter);

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
				      *ss, *mm, 1, gfn, gsbval, traps, Tsk, mask, *minprob);
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

	    /* display pmix
	       for (x=0; x<*nmix; x++) {
	       Rprintf("%3d", x);
	       for (n=0; n<*nc; n++) {                    
	       Rprintf("%5.2f", pmix[*nmix * n + x]);
	       }
	       Rprintf("\n");
	       }
	    */

            /* Loop over individuals... */
            for (n=0; n<*nc; n++) {                      /* CH numbered 0 <= n < *nc */
                a[n] = 0;
                tempsum = 0;
                if (knownclass[n] == 1) {
                    /* unknown class : weighted by probability of membership */
		    for (x=0; x<*nmix; x++) {
			pmixnx = pmixn[*nmix * n + x];
			if (pmixnx > 1e-6) {
			    temp = 0;
			    if (indiv > 0)
				asum[x] = 0;
			    for (m=0; m<*mm; m++) {
				if (indiv > 0)
				    asum[x] += pndot (m, n, 1, *ss, x, *nc, PIA0, gk0, *detect, 
		  		        *binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, *param);
				prwi = prwfn (m, n, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, 
					      detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix,
					      gfn, gsbval, traps, Tsk, mask, *minprob);
				temp += prwi;
			    }
			    a[n] += pmixnx * asum[x];
			    tempsum += pmixnx * temp;
			}
		    }    /* end of loop over mixtures */
		    templog = log(tempsum/a[n]);
		}
		else {
                    /* known class does not require probability of membership */
		    x = imax2(0, knownclass[n]-2);
		    temp = 0;
		    if (indiv > 0)
			asum[x] = 0;
		    for (m=0; m<*mm; m++) {
			if (indiv > 0)
			    asum[x] += pndot (m, n, 1, *ss, x, *nc, PIA0, gk0, *detect, 
					*binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, *param);
			prwi = prwfn (m, n, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, 
				      detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix,
				      gfn, gsbval, traps, Tsk, mask, *minprob);
			temp += prwi;
		    }
		    a[n] = asum[x];
		    templog = log(temp/a[n]);
		}
		a[n] = *area * a[n];		
                if (!R_FINITE(templog)) *resultcode = 9;
		if (*resultcode == 9) return;
		*value += templog;
                R_CheckUserInterrupt();
            }        /* end of loop over individuals */

	    /* multinomial probability of known class membership (excludes coeficient) */

	    if (known>0) {
		tempsum = 0;
		for (x=0; x<*nmix; x++) {
		    asum[x] = 0;
		    for (m=0; m<*mm; m++) {
			asum[x] += pndot (m, 0, 1, *ss, x, *nc, PIA0, gk0, *detect, 
					  *binomN, Tsk, *ss, nk, *cc0, *nmix, gsb0val, *param);
		    }
		    tempsum += asum[x] * pmixg[x];
		}
		for (x=0; x<*nmix; x++) 
		    *value += nm[x+1] * log(asum[x] * pmixg[x] / tempsum); 
		if (!R_FINITE(*value)) *resultcode = 9;
	    }

        }
    }
    /*-------------------------------------------------------------------------------------------*/

    else {  /* *like==0,2  Full or Partial likelihood */
	sumD = (double *) R_alloc(*gg, sizeof(double));
	sumDp = (double *) R_alloc(*gg, sizeof(double));
	for (g=0; g<*gg; g++) {
	    sumD[g] = 0;
	    sumDp[g] = 0;
	}
	if ((known>0) && (*nmix>1)) {
	    /* no groups, only mixture classes; use group 0 for aggregates */
	    for (m=0; m<*mm; m++)  {
		sumD[0] += Dmask[m];
		for (x=0; x<*nmix; x++) {
		    pdt = pndot (m, 0, 1, *ss, x, *gg, PIA0, gk0, *detect, *binomN, 
				 Tsk, *ss, nk, *cc0, *nmix, gsb0val, *param);
		    sumDp[0] += pdt * pmixg[x] * Dmask[m];  
		}
	    }
	}
	else {
	    /* groups incompatible with mixtures containing some known */
	    for (g=0; g<*gg; g++) {
		for (m=0; m<*mm; m++)  {
		    sumD[g] += Dmask[*mm * g + m];
		    for (x=0; x<*nmix; x++) {
			pdt = pndot (m, g, 1, *ss, x, *gg, PIA0, gk0, *detect, *binomN, 
				     Tsk, *ss, nk, *cc0, *nmix, gsb0val, *param);
			sumDp[g] += pdt * pmixg[*nmix * g + x] * Dmask[*mm * g + m];
		    }
		}
	    }
	}

	/* 3 */
	if (timing) ticks = timestamp(ticks, &counter);
                   
        *value = 0;
        /* compute likelihood component from pr(wi) */
        if (*like == 0)   /* Full likelihood only */
        for (n=0; n < *nc; n++) {
	    g = grp[n]-1;
	    temp = 0;

	    if (knownclass[n] == 1) {
		/* unknown class : weighted by probability of membership */
		for (x=0; x<*nmix; x++) {
		    pmixnx = pmixg[*nmix * g + x];
		    for (m=0; m<*mm; m++) {
			prwi = prwfn (m, n, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, 
				      detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, 
				      gfn, gsbval, traps, Tsk, mask, *minprob);
			
			tempp = prwi * pmixnx * Dmask[*mm * g + m];	    
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
	    }
	    else {
                x = imax2(0, knownclass[n]-2);
		for (m=0; m<*mm; m++) {
		    prwi = prwfn (m, n, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, 
				  detspec, h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, 
				  gfn, gsbval, traps, Tsk, mask, *minprob);		    
		    tempp = prwi * Dmask[m];	    
		    temp += tempp;
		}
	    }
	    templog = log(temp);
	    if (!R_FINITE(templog)) *resultcode = 9;
	    if (*resultcode == 9) return;
	    *value += templog;
	    R_CheckUserInterrupt();
	}

	/* 4 */
	if (timing) ticks = timestamp(ticks, &counter);

	if ((known>0) && (*nmix>1)) {
	    /* likelihood component due to n */
	    tempN = sumD[0] * *area;
            /* Poisson */
            if (*distrib == 0) *value += gpois(*nc, sumDp[0] * *area, 1);
            /* binomial */
            if (*distrib == 1) {
		if (*nc > tempN) {
		    *resultcode = 8;
		    return;
		}
		*value += gbinomFP (*nc, tempN, sumDp[0] / sumD[0], 1);
	    }
            /* likelihood component due to denominator of f() */
            *value -= *nc * log(sumDp[0]);
            /* adjustment for mixture probabilities when class known */
	    for (x=0; x<*nmix; x++)
		*value += nm[x+1] * log(pmixg[x]);
	}
	else {
	    /* no known mixture classes */
	    for (g=0; g<*gg; g++) {
		/* likelihood component due to n */
		tempN = sumD[g] * *area;
		/* Poisson */
		if (*distrib == 0) *value += gpois(ng[g], sumDp[g] * *area, 1);
		/* binomial */
		if (*distrib == 1) {
		    if (ng[g] > tempN) {
			*resultcode = 8;
			return;
		    }
		    *value += gbinomFP (ng[g], tempN, sumDp[g] / sumD[g], 1);
		}
		/* likelihood component due to denominator of f() */
		*value -= ng[g] * log(sumDp[g]);
		
		/* special case: superbinomial (specified N) */
		if (*distrib>=2) *value += gbinomFP (ng[g], (double) *distrib,
						     sumDp[g] * *area / *distrib, 1);
	    }
	}
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
  int    *detect,  /* code 0 = multicatch, 1 = proximity */
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
    int    *knownclass,  /* known membership of 'latent' classes */
    double *traps,       /* x,y locations of traps (first x, then y) */
    double *Tsk,         /* nk x s usage matrix */
    double *mask,        /* x,y points on mask (first x, then y) */

    double *gsbval,      /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) */
    int    *cc,          /* number of g0/sigma/b combinations  */
    int    *PIA,         /* lookup which g0/sigma/b combination to use for given n, S, K */

    double *area,        /* area associated with each mask point (ha) */
    double *miscparm,    /* miscellaneous parameters (cuerate etc.) */
    int    *normal,      /* code 0 don't normalise, 1 normalise */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *minprob,     /* minimum value of P(detection history) */
    double *value,       /* return values */
    int    *resultcode   /* 0 if OK */
)
{
    int    i,j,n,g,k,m,c,s,x;
    int    *ng;       /* number per group */
    int    gi;
    double *pmixg = NULL;
    double *pmixn = NULL;
    double temp;
    double sumprwi = 1.0;
    double prwi;

    double *gk;
    double *gkx;
    double *detspec;
    double *detspecx;

    int    *hc0 = NULL;        
    int    *hindex = NULL;     
    double *h = NULL;          

    int    *hc0x;
    int    *hindexx;
    double *hx;            /* 2011-11-15 */ 
    int    *start = NULL;

    gfnptr gfn;
    prwfnptr prwfn;
    int cumk[maxnpoly];
    int nk = 0;
    int nd = 0;
    int nc1 = 0;
    int nv; 
    int gpar = 2;    /* number of 'g' (detection) parameters */

    int nested = 0;  /* 0 = not nested, 1 = nested (CL only) */
    int *cumng;  /* used for 'nested' */
    double tempg;

    int    next = 0;
    int    c0 = 0;
    double p;
    double par[4];   /* passing parameter values to integr fn  */
    double stdint = 1;
    double *ex;

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

    ng = (int *) R_alloc(*gg, sizeof(int));
    fillng(*nc, *gg, grp, ng);                                    /* number per group */

    nk = fillcumk(detect, kk, cumk);                            /* detections per polygon */
    gfn = getgfn(*fn);                                          /* see utils.c */
    prwfn = getprwfn(*detect, *param);

    if ( (*detect == 3) || (*detect==4) ) 
	start = (int *) R_alloc(nc1 * *ss, sizeof(int));
    else if (((*detect>=5) && (*detect<=9)) || (*detect==12)) 
	start = (int *) R_alloc(nc1 * *ss * nk, sizeof(int));
    nd = getstart(detect, start, nc1, nc, ss, nk, w);

    pmixn = (double *)  R_alloc(nc1 * *nmix, sizeof (double));
    pmixg = (double *)  R_alloc(*gg * *nmix, sizeof (double));
    for (i=0; i < *gg * *nmix; i++) pmixg[i] = 1; /* default */
    for (i=0; i < nc1 * *nmix; i++) pmixn[i] = 1; /* default */
    gpar = getpmix(gpar, nc1, nmix, knownclass, nc, cc, ss, nk, grp, PIA, 
		   gsbval, pmixg, pmixn);

    nv = nval(detect, nc1, cc, ss, nk);
    detspec = (double *) R_alloc(nv, sizeof(double));

    gk = (double *) S_alloc(*cc * nk * *mm, sizeof(double));    /* S_alloc sets to zero */

    precompute(detect, fn, binomN, kk, mm, cc, nk, cumk,  
	       traps, mask, gsbval, miscparm, detspec, gk); 

    if (*normal && ((*detect==0) || (*detect==3) || (*detect==4))) {
        hc0 = (int *) R_alloc (*cc, sizeof(int));
        hindex = (int *) S_alloc (nc1 * *ss, sizeof(int));
	h = (double *) S_alloc (nc1 * *ss * *mm * *nmix, sizeof(double)); 
	geth (nc1, cc, nmix, nk, ss, mm, PIA, hc0, gk, Tsk, h, hindex);
    } 

    getdetspec (detect, fn, nc, nc1, cc, nmix, nd, nk, ss,     /* complete filling of detspec */
		kk, mm, PIA, miscparm, start, detspec);

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
        sumprwi = 0;
        for (x=0; x<*nmix; x++) {
            temp = 0;
            for (m=0; m<*mm; m++) {
		prwi = prwfn (m, *which-1, 1, *ss, x, w, xy, signal, PIA, gk, *binomN, detspec,
			      h, hindex, *cc, *nc, nk, *ss, *mm, *nmix, gfn, gsbval, traps, 
                              Tsk, mask, *minprob);
		temp += prwi;
            }
            sumprwi += pmixg[x] * temp;
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
    detspecx = (double *) R_alloc(nv, sizeof(double));
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
	ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
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
			   cumk[k], cumk[k+1]-1, cumk[nk], *xx, ex) / stdint;
                }
            }
        }
    }
    else if ((*detect == 4) || (*detect == 7)) {
	ex = (double *) R_alloc(10 + 2 * maxvertices, sizeof(double));
        for (c=0; c<*cc; c++) {
            par[0] = gsbval[c];
            par[1] = gsbval[*cc + c];
            par[2] = gsbval[2* *cc + c];
            stdint = gintegral1(*fn, par);
            detspecx[2+c] = stdint;              /* passed to prwitransect */
            for (k=0; k<nk; k++) {               /* over transects */
                for (m=0; m<*xx; m++) {
                    gi = i3(c,k,m,*cc,nk);
                    gkx[gi] = par[0] * integral1D (*fn, m, c, gsbval, *cc,
		       traps, X, cumk[k], cumk[k+1]-1, cumk[nk], *xx, ex) / stdint;
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
                                   p += pmixg[*nmix * g + x] * gkx[gi];
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
			      gsbval, traps, Tsk, X, *minprob);
		tempg *= prwi;
            }
            temp += tempg;
        }
        else {
            for (x=0; x<*nmix; x++) {
                temp += pmixg[x] * prwfn (i, *which-1, 1, *ss, x, w, xy, signal, PIA,
					 gkx, *binomN, detspecx, hx, hindexx, *cc, *nc, nk, 
					 *ss, *xx, *nmix, gfn, gsbval, traps, Tsk, X, *minprob);
            }
        }
        value[i] = temp / sumprwi;
    }

    *resultcode = 0;   /* successful termination pwuniform */
}
/*==============================================================================*/

