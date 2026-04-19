#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP mtgrsbed_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mtgrsbed_f(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"mtgrsbed_c", (DL_FUNC) &mtgrsbed_c, 7},
  {"mtgrsbed_f", (DL_FUNC) &mtgrsbed_f, 7},
  {NULL, NULL, 0}
};

void R_init_qgtools(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
