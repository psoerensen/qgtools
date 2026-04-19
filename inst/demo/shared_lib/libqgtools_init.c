#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call functions */
extern SEXP lib_mtgrsbed_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lib_mtgrsbed_f(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP libqgtools_version();

/* Register */
static const R_CallMethodDef CallEntries[] = {
  {"lib_mtgrsbed_c", (DL_FUNC) &lib_mtgrsbed_c, 7},
  {"lib_mtgrsbed_f", (DL_FUNC) &lib_mtgrsbed_f, 7},
  {"libqgtools_version", (DL_FUNC) &libqgtools_version, 0},
  {NULL, NULL, 0}
};

void R_init_libqgtools(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}

