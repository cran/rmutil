#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ddb(void *, void *, void *, void *, void *, void *, void *);
extern void ddp(void *, void *, void *, void *, void *, void *, void *);
extern void dmb(void *, void *, void *, void *, void *, void *, void *);
extern void dmp(void *, void *, void *, void *, void *, void *, void *);
extern void dpvfp(void *, void *, void *, void *, void *, void *, void *);
extern void inthp(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pdb(void *, void *, void *, void *, void *, void *);
extern void pdp(void *, void *, void *, void *, void *, void *);
extern void pginvgauss_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pmb(void *, void *, void *, void *, void *, void *);
extern void pmp(void *, void *, void *, void *, void *, void *);
extern void ppowexp_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ppvfp(void *, void *, void *, void *, void *, void *);
extern void psimplex_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void romberg(void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(gettvc_f)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"ddb",          (DL_FUNC) &ddb,           7},
  {"ddp",          (DL_FUNC) &ddp,           7},
  {"dmb",          (DL_FUNC) &dmb,           7},
  {"dmp",          (DL_FUNC) &dmp,           7},
  {"dpvfp",        (DL_FUNC) &dpvfp,         7},
  {"inthp",        (DL_FUNC) &inthp,         9},
  {"pdb",          (DL_FUNC) &pdb,           6},
  {"pdp",          (DL_FUNC) &pdp,           6},
  {"pginvgauss_c", (DL_FUNC) &pginvgauss_c, 10},
  {"pmb",          (DL_FUNC) &pmb,           6},
  {"pmp",          (DL_FUNC) &pmp,           6},
  {"ppowexp_c",    (DL_FUNC) &ppowexp_c,    10},
  {"ppvfp",        (DL_FUNC) &ppvfp,         6},
  {"psimplex_c",   (DL_FUNC) &psimplex_c,   10},
  {"romberg",      (DL_FUNC) &romberg,       9},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {"gettvc_f", (DL_FUNC) &F77_NAME(gettvc_f), 15},
  {NULL, NULL, 0}
};

void R_init_rmutil(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
