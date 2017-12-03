#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void alongtransectC(void *, void *, void *, void *, void *, void *, void *);
extern void chat(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fxIHP(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getdenomext(void *, void *, void *, void *, void *, void *, void *, void *);
extern void hdotpoint(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void inside(void *, void *, void *, void *, void *, void *);
extern void integralprw1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void makelookup(void *, void *, void *, void *, void *, void *, void *);
extern void naivecap2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void naived(void *, void *, void *, void *, void *, void *, void *, void *);
extern void naiveRPSV(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void nearestC(void *, void *, void *, void *, void *, void *);
extern void ontransect(void *, void *, void *, void *, void *, void *, void *);
extern void pdotpoint(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hdotpoly(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void secrloglik(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void simdetect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingsingle(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingmulti(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingproximity(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingcount(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingcapped(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingpolygon(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingpolygonX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingsignal(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingtelemetry(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingtransect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void trappingtransectX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void unmarkedloglik(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void presenceloglik(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void presenceloglik2017(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"alongtransectC",    (DL_FUNC) &alongtransectC,     7},
    {"chat",              (DL_FUNC) &chat,              29},
    {"fxIHP",             (DL_FUNC) &fxIHP,             36},
    {"getdenomext",       (DL_FUNC) &getdenomext,        8},
    {"hdotpoint",         (DL_FUNC) &hdotpoint,         13},
    {"inside",            (DL_FUNC) &inside,             6},
    {"integralprw1",      (DL_FUNC) &integralprw1,      25},
    {"makelookup",        (DL_FUNC) &makelookup,         7},
    {"naivecap2",         (DL_FUNC) &naivecap2,         11},
    {"naived",            (DL_FUNC) &naived,             8},
    {"naiveRPSV",         (DL_FUNC) &naiveRPSV,          9},
    {"nearestC",          (DL_FUNC) &nearestC,           6},
    {"ontransect",        (DL_FUNC) &ontransect,         7},
    {"pdotpoint",         (DL_FUNC) &pdotpoint,         14},
    {"hdotpoly",          (DL_FUNC) &hdotpoly,          12},
    {"secrloglik",        (DL_FUNC) &secrloglik,        42},
    {"simdetect",         (DL_FUNC) &simdetect,         28},
    {"trappingsingle",    (DL_FUNC) &trappingsingle,    18},
    {"trappingmulti",     (DL_FUNC) &trappingmulti,     18},
    {"trappingproximity", (DL_FUNC) &trappingproximity, 18},
    {"trappingcount",     (DL_FUNC) &trappingcount,     18},
    {"trappingcapped",    (DL_FUNC) &trappingcapped,    18},
    {"trappingpolygon",   (DL_FUNC) &trappingpolygon,   19},
    {"trappingpolygonX",  (DL_FUNC) &trappingpolygonX,  17},
    {"trappingsignal",    (DL_FUNC) &trappingsignal,    21},
    {"trappingtelemetry", (DL_FUNC) &trappingtelemetry, 16},
    {"trappingtransect",  (DL_FUNC) &trappingtransect,  19},
    {"trappingtransectX", (DL_FUNC) &trappingtransectX, 17},
    {"unmarkedloglik",    (DL_FUNC) &unmarkedloglik,    12},
    {"presenceloglik",    (DL_FUNC) &presenceloglik,    13},
    {"presenceloglik2017", (DL_FUNC) &presenceloglik2017, 15},
    {NULL, NULL, 0}
};

void R_init_secr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
