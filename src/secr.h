/* functions to be exported */
void gxy (
    int *n, 
    int *fn,
    double *par,
    double *w, 
    double *xy
); 

void integralprw1 (
    int    *detect,    /* detector 0 multi, 1 proximity etc. */
    double *gsb0val,   /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    int    *nc,        /* number of individuals */   
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    int    *nmix,      /* number of mixtures */ 
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *cc0,       /* number of g0/sigma/b combinations [naive animal] */
    int    *gsb0,      /* lookup which g0/sigma/b combination to use for given n, S, K [naive animal] */
    int    *ncol,      /* number of columns in gsb0; added 2009 06 25 */
    double *area,      /* area associated with each mask point (ha) */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 10 = signal strength, 11 = binary signal strength */
    int    *binomN,    /* number of trials for 'count' detector modelled with binomial */
    double *cut,       /* transformed signal strength threshold for detection */
    double *a,         /* return value integral of pr(w0) */
    int    *resultcode /* 0 for successful completion */
);


void secrloglik (
    int    *like,        /* likelihood 0 full, 1 conditional */
    int    *detect,      /* detector 0 multi, 1 proximity etc. */
    int    *distrib,     /* distribution of nc 0 Poisson, 1 binomial */
    int    *w,           /* capture histories (1:nc, 1:s, 1:k) */
    double *xy,          /* xy coordinates of polygon records */
    double *signal,      /* signal strength vector */
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
    int    *gsb,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *gsb0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive animal] */
    double *area,        /* area associated with each mask point (ha) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *cut,         /* transformed signal strength threshold for detection */
    double *minprob,     /* minimum value of P(detection history) */
    double *a,           /* a(theta) */
    double *value,       /* return value integral of pr(w0) */
    int    *resultcode   /* 0 if OK */
);

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
    int    *gsb,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *gsb0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive animal] */
    double *pID,         /* Parameter value - probability marked animal identified on resighting */
    double *area,        /* area associated with each mask point (ha) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard rate, 2 = exponential etc. */
    int    *binomN,      /*  */
    double *minprob,     /* minimum value of P(detection history) */
    double *value,       /* return value integral of pr(w0) */
    int    *resultcode   /* 0 if OK */
);

void simsecr (
    int    *detect,     /* detector 0 multi, 1 proximity */
    double *gsb0val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [naive animal] */
    double *gsb1val,    /* Parameter values (matrix nr= comb of g0,sigma,b nc=3) [caught before] */
    int    *cc0,        /* number of g0/sigma/b combinations for naive animals */
    int    *cc1,        /* number of g0/sigma/b combinations for caught before */
    int    *gsb0,       /* lookup which g0/sigma/b combination to use for given g, S, K [naive animal] */
    int    *gsb1,       /* lookup which g0/sigma/b combination to use for given n, S, K [caught before] */
    int    *N,          /* number of animals */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *nmix,       /* number of classes */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    int    *used,       /* ss x kk array of 0/1 codes for used */
    int    *Markov,     /* code 0 if behavioural response is learned, 1 if Markov */
    int    *binomN,     /* number of trials for 'count' detector modelled with binomial */
    double *cut,        /* detection threshold on transformed scale */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential, 3 = uniform */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    double *detectedXY, /* x,y locations of detections  */
    double *signal,     /* vector of signal strengths, one per detection */
    int    *value,      /* return value array of trap locations n x s */
    int    *resultcode
);

double pfn (
    int fn,
    double d2val,
    double g0,
    double sigma,
    double z,
    double area,
    double w2);

void pdotpoint (double *xy, int *nxy, double *traps, int *used, int *kk,
		int *fn, double *par, int *nocc, double *w2, double *value);

void pdotpoly (double *xy, int *nxy, double *traps, int *used, int *nk, 
    int *kk, int *fn, double *par, int *nocc, double *value);

void pdottransect (double *xy, int *nxy, double *traps, int *used, int *nk, 
    int *kk, int *fn, double *par, int *nocc, double *value);

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
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *maxone,      /* maximum of one detection per animal per occasion */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
);

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
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,          /* truncation radius */
    int    *maxone,      /* maximum of one detection per animal per occasion */
    int    *n,           /* number of individuals detected */
    int    *caught,      /* caught in session */
    double *detectedXY,  /* x,y locations of detections  */
    int    *value,       /* return value matrix of trap locations n x s */
    int    *resultcode
);

void integral2Dtest
    (int *fn, int *m, int *c, double *gsbval, int *cc, double *traps, 
    double *mask, int *n1, int *n2, int *kk, int *mm, double *result);

void trappingsignal (
    double *beta0,     /* Parameter : intercept */
    double *beta1,     /* Parameter : slope */
    double *sdS,       /* Parameter : error sd */
    double *cut,       /* detection threshold on transformed scale */
    double *sdM,       /* movement between occasions */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    int    *used,      /* ss x kk array of 0/1 codes for usage */
    int    *fn,        /* code 10 = signal strength, 11 = signal strength with spherical spr */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *signal,    /* signal strength, one per detection */  
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
);

void trappingtimes (
    double *g0,        /* Parameter : detection intercept  */
    double *sigma,     /* Parameter : detection scale */
    double *z,         /* Parameter : detection shape (hazard) */
    int    *ss,        /* number of occasions */
    int    *kk,        /* number of traps */
    int    *N,         /* number of animals */
    double *animals,   /* x,y points of animal range centres (first x, then y)  */
    double *traps,     /* x,y locations of traps (first x, then y)  */
    int    *used,      /* ss x kk array of 0/1 codes for usage */
    int    *fn,        /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,        /* truncation radius */
    int    *n,         /* number of individuals caught */
    int    *caught,    /* caught in session */
    double *times,     /* time of detection within occasion */  
    int    *value,     /* return value matrix of trap locations n x s */
    int    *resultcode
);

void trappingcount (
    double *g0,         /* Parameter : detection intercept */
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    int    *used,      /* code 0 unused 1 used for each trap x occasion */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 11 = normal signal */
    double *w2,        /* truncation radius */
    int    *binomN,     /* number of trials for 'count' detector modelled with binomial */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode
);

void trappingmulti (
    double *g0,         /* Parameter : detection magnitude */
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    int    *used,      /* code 0 unused 1 used for each trap x occasion */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,         /* truncation radius */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode  /* 0 for successful completion */
);

void trappingsingle (
    double *g0,         /* Parameter : detection magnitude */
    double *sigma,      /* Parameter : detection scale */
    double *b,          /* Parameter : detection shape (hazard) */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *N,          /* number of animals */
    double *animals,    /* x,y points of animal range centres (first x, then y)  */
    double *traps,      /* x,y locations of traps (first x, then y)  */
    int    *used ,      /* code 0 unused 1 used for each trap x occasion */
    int    *fn,         /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    double *w2,         /* truncation radius */
    int    *n,          /* number of individuals caught */
    int    *caught,     /* caught in session */
    int    *value,      /* return value matrix of trap locations n x s */
    int    *resultcode  /* 0 for successful completion */
);

void naived (
    double *sigma,      /* Parameter : detection scale */
    int    *kk,         /* number of traps */
    int    *nc,
    double *traps,      /* x,y locations of traps (first x, then y) */
    double *animals,    /* x,y locations of traps (first x, then y) */
    int    *fn,         /* code 0 = halfnormal ONLY */
    double *value       /* return value  */
);

void naiveRPSV (
    double *sigma,      /* Parameter : detection scale */
    double *z,          /* parameter : detection shape (probably fixed) */
    int    *kk,         /* number of traps */
    int    *nc,
    double *traps,      /* x,y locations of traps (first x, then y)   */
    double *animals,    /* x,y locations of animals (first x, then y) */
    int    *fn,         /* code 0 = halfnormal ONLY */
    double *value       /* return value */
);

void naivecap2 (
    int    *detect,     /* code 0 = multicatch, 1 = proximity */
    double *g0,         /* Parameter : detection magnitude */
    double *sigma,      /* Parameter : detection scale */
    int    *ss,         /* number of occasions */
    int    *kk,         /* number of traps */
    int    *mm,
    double *traps,      /* x,y locations of traps (first x, then y) */
    double *mask,       /* x,y points on mask (first x, then y) */
    int    *fn,         /* code 0 = halfnormal ONLY */
    double *value       /* return value  */
);

void makelookup (   
    double *x,            /* input matrix */
    int    *nrow,         /* input */
    int    *ncol,         /* input */
    int    *unique,       /* output number of unique rows */
    double *y,            /* output matrix of unique rows (byrow=T) */
    int    *index,        /* output lookup rows of x in y */
    int    *resultcode    /* zero if OK */
);

void nearest (   
    int    *nxy,        /* number of input points */
    double *xy,         /* input points */
    int    *ntrap,      /* input */
    double *traps,      /* input */
    int    *p,          /* output index of nearest point */
    double *d           /* output distance to nearest point */
);

void inside (double *xy, int *n1, int *n2, int *np, double *poly, int *in);

void ontransect (double *xy, int *n1, int *n2, int *npts, double *transect, double *tol, int *on);

void secrcellprob (
    int    *detect,      /* detector 0 multi, 1 proximity */
    int    *w,           /* capture histories (1:nc, 1:s, 1:k) */
    double *xy,          /* xy coordinates */
    double *signal,      /* signal strength */
    int    *grp,         /* group number for 0<=n<*nc   [full likelihood only] */
    int    *nc,          /* number of capture histories */
    int    *ss,          /* number of occasions */
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
    int    *gsb,         /* lookup which g0/sigma/b combination to use for given n, S, K */
    int    *gsb0,        /* lookup which g0/sigma/b combination to use for given g, S, K [naive animal] */
    double *area,        /* area associated with each mask point (ha) */
    int    *fn,          /* code 0 = halfnormal, 1 = hazard, 2 = exponential */
    int    *binomN,      /* number of trials for 'count' detector modelled with binomial */
    double *cut,         /* transformed signal strength threshold for detection */
    double *minprob,     
    double *value,       /* return value integral of pr(w0) */
    int    *resultcode   /* 0 if OK */
);
