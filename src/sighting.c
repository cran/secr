/* Functions related to sighting likelihood October-December 2015 */
/* tweaked 2016-10-10, revised 2016-11-18 */

/* EH eq refers Efford & Hunter ms 2016-11-18 */

/*
like 0 full likelihood
like 1 conditional likelihood
like 3 unused here
like 4 unused here
like 5 all sighting, n0 known
like 6 all sighting, n0 unknown
*/

/*
   can compile with gcc 4.6.3 :
   gcc -Ic:/R/R-3.3.0/include -c sighting.c -Wall -pedantic -std=gnu99

*/
#include "secr.h"

/*==============================================================================*/

/* set pID for each occasion and latent class */
/* called by secrloglik, fxIHP, chat */
int markresightini (
    int    ss,          /* number of occasions */
    int    nmix,        /* number of mixtures */
    int    markocc[],   /* which are marking occasions? */
    int    nk,          /* number of traps */
    int    ncol,        /* number of columns in PIA */
    int    PIA[],
    int    cc,
    double gsbval[],
    double pID[],        /* UPDATED on exit */
    int    gpar          
    ) {

    int c,k,s,x;
    int extracol = 0;

    /* set pID */
    for (s=0; s < ss; s++) {
	for (x=0; x < nmix; x++) {
	    if (markocc[s] == -1)            /* muddled sighting occasions */
		pID[s + x * ss] = 0.0;  
	    else if (markocc[s] == 0) {
		/* pick pID up from gsbval for n=0, s, k, x */
		/* column gpar+1 */
                /* find a detector that was used */
		k = 0;
		while (k<nk && PIA[i4(0, s, k, x, ncol, ss, nk)]==0) k++;
		c = PIA[i4(0, s, k, x, ncol, ss, nk)] - 1;
		/* Rprintf("c %3d gpar %3d cc %3d gsbval[c] %8.4f\n", 
		   c, gpar, cc,  gsbval[c + gpar * cc]); */
		if (c>=0)
		    pID[s + x * ss] = gsbval[c + gpar * cc]; 
		else
		    pID[s + x * ss] = 0;
		extracol = 1;
	    }
	    else pID[s + x * ss] = 1.0;  /* markocc -2, 1 - perfect ID! */  
	}
    }

    /* update gpar if pID */
    return(gpar + extracol);  /* gpar is used later if estimating pmix */
}
/*==============================================================================*/

void getpdots (int m, int n, int markocc[], int x, int ncol,
	       int PIA0[], double gk0[], double hk0[], int detect[], int binomN[], double Tsk[],
	       int ss, int nk, int cc0, int nmix, double gsb0val[], double pdots[])

/*
    probability animal at point m on mask was caught before j
    2015-10-05 adapted from pndot; revised 2016-11-19 

    See EH eq 5

    works only for like < 5 (i.e. not allsighting)
    allsighting implies caught before first occasion, so not used then

*/

{
    int k,s,c;
    int gi;
    int wxi;
    double pp;
    double p0;
    double g1;
    double Tski;

    pp = 1;
    for (s=0; s<ss; s++) {
	if (markocc[s] > 0) {  /* marking occasions only */
	    if (binomN[s] < 0) error ("negative binomN not allowed in C fn getpdots");

	    for (k=0; k< nk; k++) {
		wxi = i4(n,s,k,x,ncol,ss,nk);
		c = PIA0[wxi] - 1;
		if (c >= 0) {    /* drops unset traps */
		    gi = i3(c,k,m,cc0,nk);
		    g1 = gk0[gi];
		    Tski = Tsk[s * nk + k];
		    /* marking occasions */
		    if (fabs(Tski-1) > 1e-10) {                 /* effort <> 1.0 */
			if ((detect[s] < 8) & (detect[s] != 5))  {
			    if (binomN[s] == 0) {
				if (detect[s] == 2)
				    p0 = exp(-Tski * hazard(g1));  /* Poisson count detector  */
				else
				    p0 = exp(-Tski * g1);          /* Poisson polygon */
			    }
			    else
				p0 = pow(1-g1, Tski * binomN[s]);     /* Binomial */
			}
			else error("no effort adjustment for detector type");
		    }
		    else {
			if (binomN[s] == 0) {
			    if (detect[s] == 2)
				p0 = exp(-hazard(g1));         /* Poisson count detector  */
			    else
				p0 = exp(-g1);                 /* Poisson polygon */
			}
			else  if (binomN[s] == 1) {
                            /* expect traps to end up here */
			    p0 = 1 - g1;                           /* Bernoulli */
			}
			else {
			    p0 = pow(1-g1, binomN[s]);                /* Binomial */
			}
		    }
		    pp *= p0;
		}	
	    }
	}
	pdots[i3(s, x, m, ss, nmix)] = (1 - pp); 

    }
}
/*===============================================================*/
   

void incmusk (double Dprwi, int n, int m, int x, int PIA[], double gk[], 
	      double hk[], int binomN[], int cc, int ncol, int nk, int ss, int nmix, 
	      double gsbval[], double Tsk[], int markocc[], int firstocc, 
	      int detect[], double tmpmusk[], int debug) {
/*
  Increment expected number of sightings of _marked_ animals at each occasion s and 
  detector k.

  Function is called for each animal n at each mask point m.

  Argument w unused and removed 2019-01-05

  Intended for like < 5 (i.e. not allsighting), but prob OK for like 5 (known-n)
  if H(m).pi(m) provided instead of Dprwi

*/

    int k,s,c;
    int gi;
    int wxi;
    double H;
    double g1;
    double Tski;

    for (s=(firstocc+1); s<ss; s++) {    /* after first capture */
	if (binomN[s] < 0) error ("negative binomN not allowed in C fn incmusk");
	if (markocc[s] <= 0) {  /* sighting occasions */
	    for (k = 0; k < nk; k++) {
		H = 0;
		wxi = i4(n,s,k,x,ncol,ss,nk);
		c = PIA[wxi] - 1;
		if (c >= 0) {    /* drops unset traps */
		    gi = i3(c,k,m,cc,nk);
		    g1 = gk[gi];
		    Tski = Tsk[s * nk + k];
		    /* marking occasions */
		    if (fabs(Tski-1) > 1e-10) {                 /* effort <> 1.0 */
			if ((detect[s] < 8) & (detect[s] != 5))  {
			    if (binomN[s] == 0) {
				if (detect[s] == 2)
				    H = Tski * hazard(g1);  /* Poisson count detector  */
				else
				    H = Tski * g1;          /* Poisson polygon */
			    }
			    else
				H =  pow(1-g1, Tski * binomN[s]);     /* Binomial */
			}
			else error("no effort adjustment for detector type");
		    }
		    else {
			if (binomN[s] == 0) {
			    if (detect[s] == 2)
				H = hazard(g1);         /* Poisson count detector  */
			    else
				H = g1;                 /* Poisson polygon */
			}
			else  if (binomN[s] == 1) {
			    /* expect traps to end up here */
			    H = g1;                           /* Bernoulli */    
			}
			else {
			    H = pow(1-g1, binomN[s]);                /* Binomial */  
			}
		    }	
		    tmpmusk[s * nk + k] += H * Dprwi;
		}
	    } /* end of k loop */
	} /* end of markocc<=0 */
    } /* end of s loop */
}
/*===============================================================*/

/* void getfirstocc2(int ss, int nk, int nc, int w[], int grp[], int detect[], */
/* 		  int firstocc[], int ntelem[], int telem[], int nzero[],  */
/* 		  int allzero[]) { */
/*     /\* Find first marking occasion of each animal *\/ */
/*     /\* if not found return ss *\/ */
/*     /\* Also return count number of 'zero history' animals *\/ */
/*     /\* defined as not detected on non-telemetry occasion *\/ */

/*     int n, s, k, wi; */
/*     for (n=0; n<nc; n++) { */
/* 	firstocc[n] = ss;  /\* ss is beyond feasible first occasions *\/ */
/* 	allzero[n] = 1; */
/* 	for (s=0; s<ss; s++) { */
/* 	    for (k=0; k< nk; k++) { */
/* 		wi = i3(n, s, k, nc, ss); */
/* 		if (abs(w[wi])>0) { */
/* 		    firstocc[n] = s; */
/* 		    if (detect[s] == 13) { */
/* 			telem[n] = 1; */
/* 			ntelem[grp[n]-1] ++; */
/* 		    } */
/* 		    else { */
/* 			allzero[n] = 0; */
/* 		    } */
/* 		    break; */
/* 		} */
/* 	    } */
/* 	    if (firstocc[n] == s) break; */
/* 	} */
/* 	nzero[grp[n]-1] += allzero[n]; */

/*         /\* /\\* if caught on telemetry occasion *\\/ *\/ */
/* 	/\* wi = i3(n, ss-1, nk-1, nc, ss); *\/ */
/* 	/\* if (abs(w[wi])>0) { *\/ */
/* 	/\*     telem[n] = 1; *\/ */
/* 	/\*     ntelem[grp[n]-1] ++; *\/ */
/* 	/\* } *\/ */

/* 	/\* if (firstocc[n] > lastocc) { *\/ */
/* 	/\*     allzero[n] = 1; *\/ */
/* 	/\*     nzero[grp[n]-1] ++; *\/ */
/* 	/\* } *\/ */
/*     } */
/* } */

void getfirstocc2(int ss, int nk, int nc, int w[], int grp[], int knownclass[],
		  int detect[], int firstocc[], int ntelem[], int telem[], 
		  int nzero[], int kcnzero[],  int allzero[]) {
    /* Find first marking occasion of each animal */
    /* if not found return ss */
    /* Also return count number of 'zero history' animals */
    /* defined as not detected on non-telemetry occasion */

    int n, s, k, wi;
    for (n=0; n<nc; n++) {
	firstocc[n] = ss;  /* ss is beyond feasible first occasions */
	allzero[n] = 1;
	for (s=0; s<ss; s++) {
	    for (k=0; k< nk; k++) {
		wi = i3(n, s, k, nc, ss);
		if (abs(w[wi])>0) {
		    firstocc[n] = s;
		    if (detect[s] == 13) {
			telem[n] = 1;
			ntelem[grp[n]-1] ++;
		    }
		    else {
			allzero[n] = 0;
		    }
		    break;
		}
	    }
	    if (firstocc[n] == s) break;
	}
	nzero[grp[n]-1] += allzero[n];
	kcnzero[knownclass[n]-1] += allzero[n];
    }
}

/*===============================================================*/

/* increment each mu_sk for for one marked individual (mult=1) */
/* or for multiple all-zero histories (mult>1) */
/* DOES NOT YET ALLOW FOR MIXTURES */
void finmusk (int ss, int nk, double tmpmusk[], double musk[], double sumDprwi, 
	      double pID[], double mult, int debug)
{
    int s;
    int k;
    if (mult>0) {
	for (s = 0; s < ss; s++) {
	    for (k = 0; k < nk; k++) {
		if (debug) Rprintf("s %4d k %4d pID %6.4f tmpmusk %8.5f mult %8.5f \n",
			s,k,pID[s],tmpmusk[s*nk+k],mult);
		musk[s * nk + k] += (1-pID[s]) * tmpmusk[s * nk + k] * 
		    mult / sumDprwi;
	    }
	}
    }
}   
/*===============================================================*/

/* Log likelihood components for sightings (both unmarked and unidentified) */
 
/* Key inputs
   -- markocc  (only sighting occasions are used)
   -- pdots    pdot on successive occasions from getpdots() (like < 5)
   -- a0       mixture-specific effective area of marking in hectares (like 1,6)
   -- pi       known distribution of marked animals (Pr(animal marked is from pixel m) (like 5,6)
   -- chat     vector of overdispersion for Tu, Tm; (1,1) if not calculated
   -- Tu       sightings of unmarked animals; either scalar or matrix K x S
   -- Tm       sightings of  marked but unidentified animals; either scalar or matrix K x S

   Outputs
   -- Tulik    (quasi) log-likelihood for Tu 
   -- Tmlik    (quasi) log-likelihood for Tm 

   Quirks
   -- no allowance yet for known class membership in Tu, Tm
*/
/*==============================================================================*/

int expectedTmTu (int like, int distrib, int TmTu, int nc, int ss, int nk,
		int cc0, int nmix, double pmix[], int mm, double D[], double pi[],
		double area, int markocc[], double pdots[], 
	    	int ncol, int PIA0[], double gk0[], double hk0[], int binomN[], int detect[], 
		  double Tsk[], int nmarked[], double a0[], double pID[], double musk[]) {

    // musk is output: expected values for Tm or Tu or Tn
    // depending on TmTu switch between Tm 0, Tu 1, Tn 2
    int c,s,k,m,x, wxi;
    double mu1 = 0.0;
    double Hskx;
    double pmarked = 0;
    double Tski;
    double A, Dmarked;
    /*------------------------------------------------------------------------------*/
    /* set working variables */
    A = mm * area;

    /*-----------------------------------------------------------------------------*/
    /* loop over occasions */
    for (s=0; s < ss; s++) {
	if (markocc[s] < 1) {     /* sighting occasions only */
	    for (k=0; k < nk; k++) {
                /* add 1 to index because element 1 is reserved for nval */
		musk[s * nk + k] = 0;                          
		Tski = Tsk[s * nk + k];         /* effort; no relation to Tusk! */
		for (x=0; x<nmix; x++) {
		    mu1 = 0;
		    wxi = i4(0,s,k,x,ncol,ss,nk);   /* n = 0, cutting corners here? */
		    c = PIA0[wxi] - 1;
		    if (c >= 0) {                   /* drops unset traps */
			for (m=0; m < mm; m++) {

			    Hskx = Tski *  hk0[i3(c,k,m,cc0,nk)]; 

			    if (like < 5) {          /* like  0,2 */
				if ((markocc[s]==0) && (s>0))
				/* otherwise -1 = unresolved and pdots not needed*/
				    pmarked = pdots[i3(s-1, x, m, ss, nmix)];
				else 
				    pmarked = 0;
				if (TmTu == 0) mu1 += D[m] * Hskx * pmarked;  /* Tm */
				else  mu1 += D[m] * Hskx * (1-pmarked);   /* Tu, Tn */
			    }
			    else if (like == 5) {   /* 'sighting only' known number marked */
				Dmarked = pi[m] * nc/area;
				if (D[m] < Dmarked) {
				    // error ("negative D-Dmarked");
				    return(51);
				}
				if (TmTu == 0)  mu1 += Dmarked * Hskx;    /* Tm */    
				else mu1 += (D[m]-Dmarked) * Hskx;    /* Tu, Tn */
			    }
			    else if (like == 6) {   /* 'sighting only' unknown number marked */
				if (a0[x] > 0) {
				    Dmarked = pi[m] * nc / area * (A/a0[x]);
				    if (D[m] < Dmarked) {
					// error ("negative D-Dmarked");
					return(51);
				    }
				    if (TmTu == 0) mu1 += Dmarked * Hskx;    /* Tm */    
				    else  mu1 += (D[m] - Dmarked) * Hskx;  /* Tu, Tn */
				}
			    }
			    else error ("unknown like");
			}  /* end m loop */
		    }
		    if (TmTu == 0) {
			/* Rprintf("mu1 %8.5f pID[s] %8.5f\n", mu1, pID[s]); */
			mu1 *= (1-pID[s]);
		    }

		    musk[s * nk + k] += mu1 * pmix[x] * area;    
		}  /* end x loop */		
	    }  /* end k loop */    
	}    
    } /* end s loop */
    /*-----------------------------------------------------------------------------*/
      
    return(0); 
}
/*==============================================================================*/

/* assume musk[] already contains cumulative hazards, so general for all *like */
int Tsightinglik (int T[], int ss, int nk, int markocc[], int ncol, 
                  int detect[], double Tsk[], double musk[], int debug, 
		  double *Tlik) {

    int s,k;
    double tempmu;
      
    /* codes for count aggregation */
    /* default is no aggregation (separate T_sk each occasion and detector) */
    int TPooled;
    int TBydetector; 
   
    int TCsk= 0;
    int nused = 0;   /* number of detectors with non zero effort */
    double summu = 0;

    int *nusedk = NULL;
    double *summuk = NULL;
    int nsight = 0;     /* number of sighting occasions */
    int firstsightocc;

    *Tlik = 0; 
    if (T[0] < 0) return(0);

    /*------------------------------------------------------------------------------*/
    /* set working variables */
    firstsightocc = ss+1;
    TPooled = (T[0] == 1);      /* TPooled == 1 if there is a single summed count    */
    TBydetector = (T[0] == nk); /* TBydetector == 1 if counts are summed by detector */  
    if (TBydetector) {
	summuk = (double *)  S_alloc(nk, sizeof (double));
	nusedk = (int *)  S_alloc(nk, sizeof (int));
    }
    if (debug>1) {
	Rprintf("TPooled %4d \n", TPooled);
	Rprintf("TBydetector %4d \n", TBydetector);
    }
    /*-----------------------------------------------------------------------------*/
    /* loop over occasions */
    for (s=0; s < ss; s++) {
	if (markocc[s] < 1) {     /* sighting occasions only */
	    nsight += 1;
	    if (s < firstsightocc) firstsightocc = s;
	    for (k=0; k < nk; k++) {
		tempmu = musk[s * nk + k];
		/* Tsk is effort; no relation to TCsk! */
		nused += Tsk[s * nk + k]>0;

                /*-----------------------------------------------------------------*/
                /* compute likelihood for this cell */
		if (!TPooled && !TBydetector) {
		    /* add 1 to index because element 1 is reserved for nval */
		    TCsk = T[s * nk + k + 1];
		    if ((TCsk>0) && (tempmu<=0)) {
			// Rprintf("TCsk = %d tempmu = %8.4f\n", TCsk, tempmu);
			// Rprintf ("zero sighting probability when T number >0\n");
		       return(53);
		    }
                    /* binary (multi, proximity, presence) */
		    if ((detect[s]<2) || (detect[s] == 11)) {
			if (TCsk>1) TCsk = 1;
			if (tempmu>0)  
			    *Tlik += dbinom(TCsk, 1, 1-exp(-tempmu), 1);
		    }
                    /* count */
		    else
			*Tlik += dpois(TCsk,  tempmu, 1);
		    if (*Tlik < -1e6) {
			// Rprintf("very negative Tlik in Tsightinglik\n");
			return(54);
		    }
		}

		if (TPooled) {
		    nused += Tsk[s * nk + k]>0;
		    if (markocc[s] == 0) {
			if (detect[s] < 2)
			    summu += 1-exp(-tempmu); /* summing probabilities? */
			else
			    summu += tempmu;
		    }
		    /* else markocc < 0, marked animals not distinguished so no increment */
		}
		if (TBydetector) {
		    nusedk[k] += Tsk[s * nk + k]>0;
		    if (markocc[s] == 0) {
			if ((detect[s]<2) || (detect[s] == 11)) 
			    summuk[k] += 1-exp(-tempmu);
			else
			    summuk[k] += tempmu;
		    }
		    /* ??? else markocc < 0, marked animals not distinguished so no increment */
		}
                /*-----------------------------------------------------------------*/
	    }  /* end k loop */    
	}    
    } /* end s loop */
    /*-----------------------------------------------------------------------------*/
    
    /* Likelihood for pooled counts only */
    if (TBydetector) {
	for (k=0; k<nk; k++) {
	    /* Increment likelihood for counts by detector */
	    /* sum over s for detector k */
	    if (debug>0) 
		Rprintf("k %4d sumT %4d summu %8.3f nused %4d \n", k, T[k+1], summu, nused);
	    if (detect[firstsightocc] < 2) {  
                /* assume sighting detector same all occasions */
                /* and p constant over occasions (using arithmetic mean here) */
		if (summuk[k]>0)
//		    *Tlik += dbinom(T[k+1],  nusedk[k], summuk[k] / nk, 1);  
		    *Tlik += dbinom(T[k+1],  nusedk[k], summuk[k] / nsight, 1);  /* 2017-03-17 */
	    }
	    else
		*Tlik += dpois(T[k+1],  summuk[k], 1);
	}
    }
    else if (TPooled) {
	/* first input is the sum over s,k */
	if (debug>0) Rprintf("sumT %4d summu %8.3f nused %4d \n", T[1], summu, nused);
	if (detect[firstsightocc] < 2)   /* assume sighting detectors same on all occasions */
	    *Tlik = dbinom(T[1],  nused, summu, 1);  /* weird use of summu - to be fixed */
	else
	    *Tlik = dpois(T[1],  summu, 1);
    } 
    return(0); 
}
/*==============================================================================*/

/* estimate overdispersion of total sighting counts Tu, Tm by simulation */
/* relate to sightinglik in secr.c */
/* 2015-10-27 */

int getm(double x, int mm, double cumprob[]){
    /* inefficient search */
    /* replace with version of 'locate' from Press et al p 101 */
    int m;
    for (m=0; m<mm; m++)
	if (cumprob[m] > x) break;
    return(m);
}
/*------------------------------------------------------------------------------*/
int locate(double x, int mm, double cumprob[]){
    /* more efficient search Press et al. 1989 p 101 */
    /* assume cumprob in ascending order */
    int ju, jm, jl;
    jl = 0;
    ju = mm;
    while (ju-jl > 1) {
	jm = (ju+jl) / 2;
	if(x > cumprob[jm]) 
	    jl = jm;
	else 
	    ju = jm;
	/* Rprintf("x %8.6f cumprob[jm] %8.6f jl %4d jm %4d ju %4d \n", 
                    x, cumprob[jm], jl, jm, ju); */
    }
    return(jl);
}
/*------------------------------------------------------------------------------*/

int discreteN (double N) {
    int tn;
    tn = (int) N;
    if (N == tn) return(tn);
    else return(tn + (unif_rand() < (N-tn)));
}
/*------------------------------------------------------------------------------*/

double getpmark(int s, int m, int x, int ss, int nmix, double pdots[]) {
    if (s==0)
	return(pdots[i3(0, x, m, ss, nmix)]);
    else 
	return( pdots[i3(s, x, m, ss, nmix)] - pdots[i3(s-1, x, m, ss, nmix)]);
}
/*------------------------------------------------------------------------------*/
    
int sightingchat (int like, int detect[], int binomN[], int nc, int ss, int nk,
		  int cc0, int nmix, double pmix[], int mm, double D[], double pimask[],
		  double area, double pID[], int markocc[], int nmark,
		  int ncol, int PIA0[], double gk0[], double hk0[], double Tsk[], 
		  double a0[], int distrib, int nsim, double pdots[], double chat[]) {
    
    int    c,i, s,k,m,r,x, wxi;
    double mu1, mu2, Hskx;
    double  musk[3];
    double Tski;
    double A;
    int    *pop;
    int    *popmarked;
    int    Nmarked;
    double pmark;
    int    temppop;
    int    N;
    double *Nm;
    double *cumprob;
    double *cumprobmkd = NULL;
    double *cumprobunmkd = NULL;
    double sumNm = 0;
    double delta;
    int    np = 0;
    int    xi[3] = {0,0,0};
    double varx[3] = {0,0,0};
    double meanx[3] = {0,0,0};
    double p[3];
    double sump[3];
    double meanp[3] = {0,0,0};
    double expectedvar[3];

    double tol = 1e-8;

/* 'online' variance algorithm
   https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance 2015-10-27 
*/

/* pimask is used only for like 5,6 (sighting only) */

    /*-----------------------------------------------------*/

    if (like == 1)  /* conditional likelihood incompatible with unresolved sightings */
	for (s=0; s<ss; s++) {
	    if (markocc[s] < 0) return(2);
	}
    
    Nm = (double *) S_alloc (mm, sizeof (double));   /* expected N in each cell */
    for (m=0; m<mm; m++) {
        Nm[m] = D[m] * area;
        sumNm += Nm[m];
    }
    
    cumprob = (double *) S_alloc (mm, sizeof (double));  /* cell membership */

    if (distrib) {  /* only for fixed N multinomial */
	cumprob[0] = Nm[0] / sumNm;
	for (m=1; m<mm; m++) {
	    cumprob[m] = cumprob[m-1] + Nm[m] / sumNm;
	}
    }

    /* sighting-only, conditioning on Nm */
    if (like==5 || like==6) {
	cumprobmkd = (double *) S_alloc (mm, sizeof (double));    /* marked cell membership */
	cumprobunmkd = (double *) S_alloc (mm, sizeof (double));  /* unmarked cell membership */
	cumprobmkd[0] = pimask[0];
	cumprobunmkd[0] = (Nm[0] - nmark * pimask[0]) / (sumNm - nmark);
	for (m=1; m<mm; m++) {
	    cumprobmkd[m] = cumprobmkd[m-1] + pimask[m];
	    cumprobunmkd[m] = cumprobunmkd[m-1] + (Nm[m] - nmark * pimask[m]) / (sumNm - nmark) ;
	}
    }

    /* for(m=0; m<mm; m+=10) Rprintf("m %5d Nm[m] %8.6f sumNm %9.4f cumprob[m] %8.6f \n", 
       m, Nm[m], sumNm, cumprob[m]); */

    pop       = (int *) R_alloc (mm*nmix, sizeof (int));  /* unmarked population in each cell */
    popmarked = (int *) R_alloc (mm*nmix, sizeof (int));  /* marked population in each cell */
    GetRNGstate();                             /* initialise random seed */
    A = mm * area;                             /* total mask area */ 

    /*-----------------------------------------------------*/

    for (r=0; r < nsim; r++) {
	for (i=0; i<3; i++) {
	    xi[i] = 0; 
	    sump[i] = 0;
	}
	np = 0;
        for(m=0; m<mm; m++) {pop[m] = 0; popmarked[m] = 0; }

        /* pre-set marked population if conditioning on nmark */
	if (like==5 || like==6) {
	    for (i=0; i<nmark; i++) {
		m = locate(unif_rand(), mm, cumprobmkd);
		popmarked[m] += 1;
	    }
	}

        /* Random initial population */

	if (distrib) {   /* N 'fixed' (= binomial n);  multinomial cells */
            /* simulated density is matched by discreteN only in the long-run average */
	    N = discreteN(sumNm);  
	    if (like==5 || like==6) {
		if (N>nmark) {
		    for (i=0; i<(N-nmark); i++) {
			m = locate(unif_rand(), mm, cumprobunmkd);
			pop[m] += 1;    /* unmarked animals in cell m */
		    }
		}
	    }
	    else {
		for (i=0; i<N; i++) {
		    m = locate(unif_rand(), mm, cumprob);
		    if ((m<0) || (m>=mm)) {
			Rprintf("erroneous location of simulated animal in sightingchat\n");
			PutRNGstate();  /* return random seed to R */
			return(1);
		    }
		    pop[m] += 1;
		}
	    }
	    /* for(m=0; m<mm; m+=10) Rprintf("sumNm %10.6f  N  %5d m %5d pop[m] %5d \n", 
	       sumNm, N, m, pop[m]); */
	}
	else {   /* Poisson total, cells multinomial allowing for marked */
	    if (like==5 || like==6) {
		N = rpois(sumNm);
		if (N>nmark) {
		    for (i=0; i<(N-nmark); i++) {
			m = locate(unif_rand(), mm, cumprobunmkd);
			pop[m] += 1;
		    }
		}
	    }
	    else {     /* Poisson total, Poisson each cell */
		for (m=0; m < mm; m++) {
		    pop[m] = rpois(Nm[m]);
		}
	    }
	}
        /* initial pop same for all mixture classes; copy if required */
	for (x=1; x<nmix; x++) {
	    pop[x*mm+m] = pop[m];
	}

	mu1 = 0; mu2 = 0; 
	for (s=0; s < ss; s++) {

	    if (markocc[s]==1) {          /* marking occasions */
		/* update marked and unmarked populations */
		for (x=0; x<nmix; x++) {
		    Nmarked=0;                   // for check
		    for (m=0; m < mm; m++) {
			temppop = pop[x*mm + m];
			if (temppop > 0) {
			    pmark = getpmark(s, m, x, ss, nmix, pdots); // difference , not cumul
			    for (i=0; i<temppop; i++) {
				if (unif_rand() < pmark) {
				    pop[x * mm + m] --;
				    popmarked[x * mm + m] ++;
				    Nmarked++;   // for check
				}
			    }
			}
		    }
		}
	    } 
	    else                          /* sighting occasions */
                /* accumulate sightings */
		for (k=0; k < nk; k++) {
		    for (i=0; i<3; i++) musk[i] = 0;
		    for (x=0; x<nmix; x++) {
			mu1 = 0; mu2 = 0; 
			wxi = i4(0,s,k,x,ncol,ss,nk);    /* n = 0, generic individual */
			c = PIA0[wxi] - 1;
			Tski = Tsk[s * nk + k];          /* effort (usage) adjustment */
			if (c >= 0) {                    /* drops unset traps */
			    for (m=0; m < mm; m++) {
				if (pop[x*mm+m]>0) {     /* any unmarked animals at m */
				    /* 2015-12-21 if Poisson convert to lambda */
				    if ((binomN[s] == 0) & (detect[s] == 2))
					Hskx = Tski * hk0[i3(c,k,m,cc0,nk)];	 
				    else
//					Hskx = Tski * gk0[i3(c,k,m,cc0,nk)];	
					Hskx = Tski * hk0[i3(c,k,m,cc0,nk)];	
	 
				    if (like == 1) {        /* CL */
					if (markocc[s] == 0) 
					    mu2 +=  popmarked[x*mm+m] * Hskx;
				    }
				    else if (like == 5) {   /* like 5 all sighting, known */
					mu1 += pop[x*mm+m] * Hskx;
					if (markocc[s] == -1)
					    mu1 += popmarked[x*mm+m] * Hskx;
					if (markocc[s] == 0) 
					    mu2 +=  pimask[m] * nc * Hskx; 
				    }
				    else if (like == 6) {   /* like 6 all sighting, unknown */
					mu1 += pop[x*mm+m] * Hskx;
					if (markocc[s] == -1)
					    mu1 += popmarked[x*mm+m] * Hskx;
					if (markocc[s] == 0) 
					    if (a0[x] > 0)   // 2015-12-31
						mu2 +=  pimask[m] * nc * (A/a0[x]) * Hskx; 
				    }
				    else {                  /* like  0,2,3,4 */
					mu1 += pop[x*mm+m] * Hskx ;
					if (markocc[s] == -1)
					    mu1 += popmarked[x*mm+m] * Hskx;
					if (markocc[s] == 0) {
					    mu2 += popmarked[x*mm+m] * Hskx;
					}
				    }
				}
			    }
			}
			musk[0] += mu1 * pmix[x];
			musk[1] += mu2 * (1 - pID[s + ss*x]) * pmix[x];
			musk[2] = musk[0] + musk[1];
		    
		    }  /* end loop over latent classes */

		    if (detect[s] < 2) {   /* multi, proximity */
			np += 1;
			for (i=0; i<3; i++) {
			    p[i] = 1-exp(-musk[i]);
			    sump[i] += p[i];
			    if (musk[i]>tol) xi[i] += (unif_rand() < p[i]);
			}
		    }
		    else {                /* count etc. */
			for (i=0; i<3; i++) {
			    if (musk[i]>tol) xi[i] += rpois(musk[i]);
			}
		    }
		}  /* end loop over detectors */    
	}  /* end loop over occasions */

	for (i=0; i<3; i++) {
	    delta = xi[i] - meanx[i];
	    meanx[i] += delta/(r+1);
	    varx[i] += delta*(xi[i] - meanx[i]);
	    /* assuming uniform probability across cells ? */
	    if (np>0) {
		delta = sump[i]/np - meanp[i];
		meanp[i] += delta/(r+1);
	    }
	}
    }  /* end loop over replicates */
    /*-----------------------------------------------------*/

    for (i=0; i<3; i++) {
	varx[i] /= nsim-1;
	if (np>0)   /* multi, proximity; implicitly assume all occasions same */
	    expectedvar[i] = np * meanp[i] * (1-meanp[i]);    
	else        /* Poisson */
	    expectedvar[i] = meanx[i];
        if (expectedvar[i] > 0) 
	    chat[i] = varx[i]/expectedvar[i]; 
	else chat[i] = 1;
    }

/*
  Rprintf("varx1 %8.4f varx2 %8.4f meanx1 %8.4f meanx2 %8.4f chat1 %8.4f chat2 %8.4f\n",
  varx1, varx2,meanx1, meanx2, chat[0], chat[1]);
*/

    PutRNGstate();  /* return random seed to R */
    return(0);
}
/*==============================================================================*/

