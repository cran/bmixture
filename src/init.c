#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bmix_gamma_fixed_k(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bmix_gamma_unknown_k(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bmix_t_fixed_k(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bmix_t_unknown_k(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dmixgamma_hat_x_seq_fixed_k(void *, void *, void *, void *, void *, void *, void *, void *);
extern void dmixgamma_hat_x_seq_unknow_k(void *, void *, void *, void *, void *, void *, void *, void *);
extern void dmixt_hat_x_seq_fixed_k(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void dmixt_hat_x_seq_unknow_k(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bmix_gamma_fixed_k",           (DL_FUNC) &bmix_gamma_fixed_k,           15},
    {"bmix_gamma_unknown_k",         (DL_FUNC) &bmix_gamma_unknown_k,         19},
    {"bmix_t_fixed_k",               (DL_FUNC) &bmix_t_fixed_k,               18},
    {"bmix_t_unknown_k",             (DL_FUNC) &bmix_t_unknown_k,             23},
    {"dmixgamma_hat_x_seq_fixed_k",  (DL_FUNC) &dmixgamma_hat_x_seq_fixed_k,   8},
    {"dmixgamma_hat_x_seq_unknow_k", (DL_FUNC) &dmixgamma_hat_x_seq_unknow_k,  8},
    {"dmixt_hat_x_seq_fixed_k",      (DL_FUNC) &dmixt_hat_x_seq_fixed_k,       9},
    {"dmixt_hat_x_seq_unknow_k",     (DL_FUNC) &dmixt_hat_x_seq_unknow_k,      9},
    {NULL, NULL, 0}
};

void R_init_bmixture(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
