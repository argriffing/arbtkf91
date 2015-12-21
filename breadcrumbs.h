#ifndef BREADCRUMBS_H
#define BREADCRUMBS_H

/*
 * Breadcrumbs for the traceback stage of dynamic programming may indicate
 * that multiple paths have identical scores.
 * To show that multiple traceback continuations from a partial solution
 * are equally valid, use a bitwise combination of these flags.
 *
 * If multiple directions are equally valid and exact equality is
 * detected symbolically, then the numerical ambiguity flag is not used.
 * The numerical ambiguity flag is used only if numerical (not symbolic)
 * ambiguity is detected.
 */
#define CRUMB_TOP  0x01
#define CRUMB_DIAG 0x02
#define CRUMB_LEFT 0x04
#define CRUMB_NUMERICAL_AMBIGUITY 0x08

typedef unsigned char breadcrumb_t;
typedef breadcrumb_t * breadcrumb_ptr;

typedef struct
{
    breadcrumb_t *data;
    slong nrows;
    slong ncols;
} breadcrumb_mat_struct;
typedef breadcrumb_mat_struct breadcrumb_mat_t[1];


#ifdef __cplusplus
extern "C" {
#endif



void breadcrumb_mat_init(breadcrumb_mat_t mat, slong nrows, slong ncols);
void breadcrumb_mat_clear(breadcrumb_mat_t mat);
breadcrumb_ptr breadcrumb_mat_entry(breadcrumb_mat_t mat, slong i, slong j);
void breadcrumb_mat_get_alignment(char **psa, char **psb,
        breadcrumb_mat_t mat, const slong *A, const slong *B);



#ifdef __cplusplus
}
#endif

#endif
