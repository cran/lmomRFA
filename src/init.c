
// Declarations for registration of native routines

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(qkap)(void *, void *, void *);
extern void F77_NAME(regtst)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"qkap",   (DL_FUNC) &F77_NAME(qkap),    3},
    {"regtst", (DL_FUNC) &F77_NAME(regtst), 17},
    {NULL, NULL, 0}
};

void R_init_lmomRFA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

