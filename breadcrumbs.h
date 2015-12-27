#ifndef BREADCRUMBS_H
#define BREADCRUMBS_H

/*
 * Breadcrumbs for the traceback stage of dynamic programming
 */

/* traceback directions for the actual alignment */
#define CRUMB_TOP  0x01
#define CRUMB_DIAG 0x02
#define CRUMB_LEFT 0x04

/* traceback directions for information propagation */
#define CRUMB_DIAG2 0x08
#define CRUMB_LEFT2 0x10

/* indicates whether max{m1, m2} or max{m0, m1, m2} or both
 * are requested to be computed for the cell */
#define CRUMB_WANT2 0x20
#define CRUMB_WANT3 0x40
#define CRUMB_CONTENDER 0x80

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
void breadcrumb_mat_get_alignment(char **psa, char **psb,
        breadcrumb_mat_t mat, const slong *A, const slong *B);
void breadcrumb_mat_get_mask(breadcrumb_mat_t mask, const breadcrumb_mat_t mat);
void breadcrumb_mat_fprint(FILE *stream, const breadcrumb_mat_t mat);
slong breadcrumb_mat_nnz(const breadcrumb_mat_t mat, breadcrumb_t pattern);

static __inline__ breadcrumb_ptr
breadcrumb_mat_srcentry(const breadcrumb_mat_t mat, slong i, slong j)
{
    return mat->data + i * mat->ncols + j;
}

static __inline__ breadcrumb_ptr
breadcrumb_mat_entry(breadcrumb_mat_t mat, slong i, slong j)
{
    return mat->data + i * mat->ncols + j;
}

static __inline__ slong
breadcrumb_mat_nrows(const breadcrumb_mat_t mat)
{
    return mat->nrows;
}

static __inline__ slong
breadcrumb_mat_ncols(const breadcrumb_mat_t mat)
{
    return mat->ncols;
}



#ifdef __cplusplus
}
#endif

#endif
