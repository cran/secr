/* Functions related to sighting likelihood October-December 2015 */

/*
   can compile with gcc 4.6.3 :
   gcc -Ic:/R/R-3.2.0/include -c sighting.c -Wall -pedantic -std=gnu99

*/
#include "secr.h"

/*==============================================================================*/

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

    int c,s,x;
    int extracol = 0;

    /* set pID */
    for (s=0; s < ss; s++) {
	for (x=0; x < nmix; x++) {
	    if (markocc[s]<0)            /* muddled sighting occasions */
		pID[s + x * ss] = 0.0;  
	    else if (markocc[s]==1)
		pID[s + x * ss] = 1.0;  /* perfect ID in the hand! */  
	    else {
		/* pick pID up from gsbval for n=0, s, k=0, x */
		/* column gpar+1 */
		c = PIA[i4(0, s, 0, x, ncol, ss, nk)] - 1;
		pID[s + x * ss] = gsbval[c + gpar * cc]; 
		extracol = 1;
	    }
	}
    }

    /* update gpar if pID */
    return(gpar + extracol);  /* gpar is used later if estimating pmix */
}
/*==============================================================================*/

void getpdots (int m, int n, int markocc[], int x, int ncol,
             int PIA0[], double gk0[], int detect, int binomN, double Tsk[], int ss,
             int nk, int cc0, int nmix, double gsb0val[], double pdots[])

/*
    probability animal at point m on mask is caught, for each critical occasion
    2015-10-05 adapted from pndot 

    works only for like < 5 (i.e. not allsighting)

*/

{
    int k,s,c;
    int gi;
    int wxi;
    double pp;
    double p0;
    double g1;
    double Tski;
    if (binomN < 0) error ("negative binomN not allowed in C fn getpdots");
	pp = 1;
	for (s=0; s<ss; s++) {
	    if (markocc[s] > 0) {
		for (k=0; k< nk; k++) {
		    wxi = i4(n,s,k,x,ncol,ss,nk);
		    c = PIA0[wxi] - 1;
		    if (c >= 0) {    /* drops unset traps */
			gi = i3(c,k,m,cc0,nk);
			g1 = gk0[gi];
			Tski = Tsk[s * nk + k];
			/* marking occasions */
			/* expect binomN = 1 if not count detector */
			if (fabs(Tski-1) > 1e-10) {                 /* effort <> 1.0 */
			    if ((detect < 8) & (detect != 5))  {
				if (binomN == 0) {
				    if (detect == 2)
					p0 = exp(-Tski * hazard(g1));  /* Poisson count detector  */
				    else
					p0 = exp(-Tski * g1);          /* Poisson polygon */
				}
				else
				    p0 = pow(1-g1, Tski * binomN);     /* Binomial */
			    }
			    else error("no effort adjustment for detector type");
			}
			else {
			    if (binomN == 0) {
				if (detect == 2)
				    p0 = exp(-hazard(g1));         /* Poisson count detector  */
				else
				    p0 = exp(-g1);                 /* Poisson polygon */
			    }
			    else  if (binomN == 1)
				p0 = 1 - g1;                           /* Bernoulli */
			    else {
				p0 = pow(1-g1, binomN);                /* Binomial */
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

int sightinglik (int Tu[], int Tm[], int like, int nc, int ss, int nk,
		  int cc0, int nmix, double pmix[], int mm, double D[], double pi[],
		  double area, double pID[], int markocc[], double pdots[], 
		 int ncol, int PIA0[], double gk0[], int binomN, int detect, 
                 double Tsk[], double a0[], double chat[], double *Tulik, double *Tmlik) {
    int c,s,k,m,x, wxi;
    double mu1, mu2, musk1, musk2, tauskx, pdot;
    double  pmarked = 0;
    double Tski, g1;
    double A, Dmarked;
    
    /* accumulators used when data not already pooled */
    int sumTu = 0;
    int sumTm = 0;
    double summu1 = 0;
    double summu2 = 0;
    
    /* for traditional measures of overdispersion */
    int nsk = 0;
    double Cp = 0;
    double S = 0;
    
    int TuPooled;    /* code to distinguish single summed count from separate Tu_sk */
    int TmPooled;    /* code to distinguish single summed count from separate Tm_sk */
    int Tusk, Tmsk;

    *Tulik = 0;    
    *Tmlik = 0; 
    
    /*------------------------------------------------------------------------------*/
    /* check inputs 2015-11-16 */
    if ((chat[0]<=0) || (chat[1]<=0)) {
	Rprintf("chat must be positive in sightinglik\n");
	return(11);
    }
    for (x=0; x<nmix; x++) {
	if ((like < 5)) {    // otherwise pdots not used
	    for (m=0; m<mm; m++)
		if (pdots[i3(ss-1, x, m, ss, nmix)] <= 0) {
		    /* Rprintf("pdot must be positive in sightinglik\n"); */
		    return(21);
		}
	}        
    }
    /*------------------------------------------------------------------------------*/
    /* set working variables */
    TuPooled = (Tu[0] == 1);   
    TmPooled = (Tm[0] == 1);   
    A = mm * area;
    /*-----------------------------------------------------------------------------*/
    /* loop over occasions */
    for (s=0; s < ss; s++) {
	// Rprintf("s %4d pID %8.6f\n", s, pID[s]); 
	if (markocc[s] < 1) {     /* sighting occasions only */
	    for (k=0; k < nk; k++) {
                /* add 1 to index because element 1 is reserved for nval */
		if (!TuPooled) {
		    Tusk = Tu[s * nk + k + 1];
		    sumTu += Tusk;
		}
		if (!TmPooled) {
		    Tmsk = Tm[s * nk + k + 1];
		    sumTm += Tmsk;
		}
		musk1 = 0; musk2 = 0;
		for (x=0; x<nmix; x++) {
		    mu1 = 0; mu2 = 0;
		    wxi = i4(0,s,k,x,ncol,ss,nk);   /* n = 0, cutting corners here? */
		    c = PIA0[wxi] - 1;
		    Tski = Tsk[s * nk + k];         /* effort; no relation to Tmsk, Tusk! */
		    if (c >= 0) {                   /* drops unset traps */
			nsk += (x == 0);            /* accumulate number of valid counts */
			for (m=0; m < mm; m++) {
			    g1 = gk0[i3(c,k,m,cc0,nk)];	
			    /* 2015-12-21 if Poisson convert to lambda */
			    if ((binomN == 0) & (detect == 2))
				tauskx = Tski * hazard(g1);   /* Poisson count detector */
			    else 
				tauskx = Tski * g1;	      /* Bernoulli, and cum hazard poly */

			    if (like == 1) {        /* CL; no Tu contribution */
				if (markocc[s] == 0)   /* otherwise -1 no marked nonID */
				    if (a0[x] > 0)   // 2015-12-31
				    mu2 +=  nc / a0[x] * tauskx * pdots[i3(s, x, m, ss, nmix)];
			    }
			    else if (like == 5) {   /* 'sighting only' known number marked */
				Dmarked = pi[m] * nc/area;
				if (D[m] < Dmarked) {
				    *Tulik = -1e20;
				    return(51);
				}
				mu1 += (D[m] - Dmarked) * tauskx;
				mu2 += Dmarked * tauskx; 
			    }
			    else if (like == 6) {   /* 'sighting only' unknown number marked */
				if (a0[x] > 0) {  // 2015-12-31
				    Dmarked = pi[m] * nc / area * (A/a0[x]);
				    if (D[m] < Dmarked) {
					*Tulik = -1e20;
					return(61);
				    }
				    mu1 += (D[m] - Dmarked) * tauskx;
				    mu2 += Dmarked * tauskx; 
				}
			    }
			    else {                  /* like  0,2,3,4 */
				if (markocc[s] == 0)   /* otherwise -1 no marked nonID */
				    pmarked = pdots[i3(s, x, m, ss, nmix)];
				else 
				    pmarked = 0;
				mu1 += D[m] * tauskx * (1 - pmarked);

                                /* divide by pdot because we know animals were marked sometime */
				if (markocc[s] == 0) {   /* otherwise -1 no marked nonID */
				    pdot = pdots[i3(ss-1, x, m, ss, nmix)]; /* ss-1 : final value */
				    mu2 += D[m] * tauskx * pmarked / pdot; 
				}
			    }
			}  /* end m loop */
		    }
                    /* mu1, mu2 both densities, so multiply by area of pixel */
		    musk1 += mu1 * pmix[x] * area;
		    musk2 += mu2 * (1-pID[s + ss*x]) * pmix[x] * area;
		    /* Rprintf("s %4d Tusk %4d musk1 %6.4f A %6.2f a0[0] %8.4f \n",
                                s,Tusk, musk1, A, a0[0]); */
		}  /* end x loop */		
		
                /*-----------------------------------------------------------------*/
                /* either increment summed mu, or compute likelihood for this cell */
		if (TuPooled) {
		    summu1 += musk1;
		}
		else if (Tu[0] >= 0) {
		    /* if ((Tusk>0) && (mu1<=0))
		       error ("zero sighting probability when number >0"); */
		    *Tulik += dpois(Tusk, musk1, 1) / chat[0];	
		    // Rprintf("s %4d Tusk %4d musk1 %6.4f Tulik %7.4f \n",
		    //	s,Tusk, musk1, *Tulik);
		}
		if (TmPooled) {
		    if (markocc[s] == 0) {
			summu2 += musk2;
		    }
		    /* else maskocc < 0, marked animals not distinguished so no increment */
		}
		else if (Tm[0] >= 0) {
		    /* if ((Tmsk>0) && (mu2<=0))
		       error ("zero sighting probability when number >0"); */
		    *Tmlik += dpois(Tmsk,  musk2, 1) /  chat[1];
		}
                /*-----------------------------------------------------------------*/
		
                /* extra code for trad measures of overdispersion */
		Cp += (Tusk - mu1) * (Tusk - mu1) / mu1;
		S  += (Tusk - mu1) / mu1;
		
	    }  /* end k loop */    
	}    
    } /* end s loop */
    /*-----------------------------------------------------------------------------*/
    
    /* ad hoc test of overdispersion measures */
    /* assuming 3 parameters estimated and at least 4 counts*/
    /*
      Cp = Cp / (nsk - 3);
      S = S / nsk; 
      Rprintf ("Chatp = %8.5f  ChatF %8.5f \n", Cp, Cp/(1+S));
    */
    
    /* Likelihood for pooled counts only */
    if ((Tu[0]>=0) && TuPooled) {
	sumTu = Tu[1];  /* first input is the sum over s,k */
	*Tulik = dpois(sumTu, summu1, 1) / chat[0];	
    }
    
    if ((Tm[0]>=0) && TmPooled) {
	sumTm = Tm[1];  /* first input is the sum over s,k */
	*Tmlik = dpois(sumTm,  summu2, 1) /  chat[1];
    }
    
    
    // Rprintf("sumTu %6d sumTm %6d summu1 %12.8f summu2 %12.8f Tulik %12.8f Tmlik %12.8f \n", 
    //   sumTu, sumTm, summu1, summu2, *Tulik, *Tmlik); 
    
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
    
int sightingchat (int like, int detect, int binomN, int nc, int ss, int nk,
		   int cc0, int nmix, double pmix[], int mm, double D[], double pi[],
		   double area, double pID[], int markocc[], 
		   int ncol, int PIA0[], double gk0[], double Tsk[], 
		  double a0[], int distrib, int nsim, double pdots[], double chat[]) {
    
    int    c,i, s,k,m,r,x, wxi;
    double mu1, mu2, musk1, musk2, tauskx;
    double Tski, g1;
    double A;
    int    *pop;
    int    *popmarked;
    int    Nmarked;
    double pmark;
    int temppop;
    int    N;
    double *Nm;
    double *cumprob;
    double sumNm = 0;
    int    x1 = 0;
    int    x2 = 0;
    double delta;
    double varx1 = 0;
    double varx2 = 0;
    double meanx1 = 0;
    double meanx2 = 0;
    double tol = 1e-8;

/* 'online' variance algorithm
   https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance 2015-10-27 
*/
    /*-----------------------------------------------------*/

    if (like == 1)  /* conditional likelihood incompatible with unresolved sightings */
    for (s=0; s<ss; s++)
	if (markocc[s] == -1) 
	    return(2);

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

    /* for(m=0; m<mm; m+=10) Rprintf("m %5d Nm[m] %8.6f sumNm %9.4f cumprob[m] %8.6f \n", 
       m, Nm[m], sumNm, cumprob[m]); */

    pop       = (int *) R_alloc (mm*nmix, sizeof (int));  /* unmarked population in each cell */
    popmarked = (int *) R_alloc (mm*nmix, sizeof (int));  /* marked population in each cell */
    GetRNGstate();                             /* initialise random seed */
    A = mm * area;                             /* total mask area */ 

    /*-----------------------------------------------------*/

    for (r=0; r < nsim; r++) {
	x1 = 0;
	x2 = 0;

        /* Random initial population */
        for(m=0; m<mm; m++) {pop[m] = 0; popmarked[m] = 0; }
	if (distrib) {   /* N 'fixed' (= binomial n);  multinomial cells */
            /* simulated density is matched by discreteN only in the long-run average */
	    N = discreteN(sumNm);  
	    for (i=0; i<N; i++) {
		m = locate(unif_rand(), mm, cumprob);
		if ((m<0) || (m>=mm)) {
		    Rprintf("erroneous location of simulated animal in sightingchat\n");
		    PutRNGstate();  /* return random seed to R */
		    return(1);
		}
		pop[m] += 1;
	    }
	    /* for(m=0; m<mm; m+=10) Rprintf("sumNm %10.6f  N  %5d m %5d pop[m] %5d \n", 
	       sumNm, N, m, pop[m]); */
	}
	else {      /* Poisson total, Poisson each cell */
	    for (m=0; m < mm; m++) {
		pop[m] = rpois(Nm[m]);
	    }
	}
        /* initial pop same for all mixture classes; copy if required */
	for (x=1; x<nmix; x++)
	    pop[x*mm+m] = pop[m];

        // For like = 5,6 2015-12-04
        // could also pre-assign population to unmarked (pop) and marked (popmarked)
        // but this would require resolution of the different distributions
        // D(x) and pi(x)

        /* Simulate detections */
	mu1 = 0; mu2 = 0; 
	for (s=0; s < ss; s++) {

            /* update unmarked population */
	    if (markocc[s]==1) {          /* marking occasions */
		for (x=0; x<nmix; x++) {
		    Nmarked=0;                   // for check
		    for (m=0; m < mm; m++) {
			if (pop[x*mm + m]>0) {
			    pmark = getpmark(s, m, x, ss, nmix, pdots);
			    temppop = pop[x*mm + m];
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

	    for (k=0; k < nk; k++) {
		musk1 = 0; musk2 = 0;
		for (x=0; x<nmix; x++) {
		    mu1 = 0; mu2 = 0; 
		    wxi = i4(0,s,k,x,ncol,ss,nk);       /* n = 0, generic individual */
		    c = PIA0[wxi] - 1;
		    Tski = Tsk[s * nk + k];             /* effort (usage) adjustment */
		    if (c >= 0) {                       /* drops unset traps */
			for (m=0; m < mm; m++) {
			    if (pop[x*mm+m]>0) {
				g1 = gk0[i3(c,k,m,cc0,nk)];	
                                /* 2015-12-21 if Poisson convert to lambda */
				if ((binomN == 0) & (detect == 2))
				    tauskx = Tski * hazard(g1);	 
				else
				    tauskx = Tski * g1;	 
				if (markocc[s] < 1) {                     /* sighting occasions */
				    if (like == 1) {        /* CL */
					if (markocc[s] == 0) 
					    mu2 +=  popmarked[x*mm+m] * tauskx;
				    }
				    else if (like == 5) {   /* like 5 all sighting, known */
					mu1 += pop[x*mm+m] * tauskx;
					if (markocc[s] == -1)
					    mu1 += popmarked[x*mm+m] * tauskx;
					if (markocc[s] == 0) 
					    mu2 +=  pi[m] * nc * tauskx; 
				    }
				    else if (like == 6) {   /* like 6 all sighting, unknown */
					mu1 += pop[x*mm+m] * tauskx;
					if (markocc[s] == -1)
					    mu1 += popmarked[x*mm+m] * tauskx;
					if (markocc[s] == 0) 
					    if (a0[x] > 0)   // 2015-12-31
						mu2 +=  pi[m] * nc * (A/a0[x]) * tauskx; 
				    }
				    else {                  /* like  0,2,3,4 */
					mu1 += pop[x*mm+m] * tauskx ;
					if (markocc[s] == -1)
					    mu1 += popmarked[x*mm+m] * tauskx;
					if (markocc[s] == 0) {
					    mu2 += popmarked[x*mm+m] * tauskx;
					}
				    }
				}
			    }
			}
		    }
		    musk1 += mu1 * pmix[x];
		    musk2 += mu2 * (1 - pID[s + ss*x]) * pmix[x];
		    
		}  /* end of x loop */
		if (musk1>tol) x1 += rpois(musk1);
		if (musk2>tol) x2 += rpois(musk2);
	    }    
	}
        delta = x1 - meanx1;
        meanx1 += delta/(r+1);
        varx1 += delta*(x1 - meanx1);

        delta = x2 - meanx2;
        meanx2 += delta/(r+1);
        varx2 += delta*(x2 - meanx2);
    }
    /*-----------------------------------------------------*/

    varx1 /= nsim-1;
    varx2 /= nsim-1;
    if (meanx1 > 0) chat[0] = varx1/meanx1; else chat[0] = 1;
    if (meanx2 > 0) chat[1] = varx2/meanx2; else chat[1] = 1;

/*
    Rprintf("varx1 %8.4f varx2 %8.4f meanx1 %8.4f meanx2 %8.4f chat1 %8.4f chat2 %8.4f\n",
       varx1, varx2,meanx1, meanx2, chat[0], chat[1]);
*/

   PutRNGstate();  /* return random seed to R */
   return(0);
}
 /*==============================================================================*/

