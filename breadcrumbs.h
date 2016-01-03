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
typedef breadcrumb_mat_struct * breadcrumb_mat_ptr;


#ifdef __cplusplus
extern "C" {
#endif


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


void breadcrumb_mat_init(breadcrumb_mat_t mat, slong nrows, slong ncols);
void breadcrumb_mat_clear(breadcrumb_mat_t mat);
void breadcrumb_mat_set(breadcrumb_mat_t mat, const breadcrumb_mat_t src);
void breadcrumb_mat_get_alignment(char *sa, char *sb, slong *plen,
        breadcrumb_mat_t mat, const slong *A, const slong *B);
void breadcrumb_mat_get_mask(breadcrumb_mat_t mask, const breadcrumb_mat_t mat);
void breadcrumb_mat_fprint(FILE *stream, const breadcrumb_mat_t mat);
void breadcrumb_mat_check_alignment(
        int *p_is_optimal, int *p_is_canonical,
        const breadcrumb_mat_t mat, const slong *A, const slong *B, slong len);

static __inline__ breadcrumb_ptr
breadcrumb_mat_srcentry(const breadcrumb_mat_t mat, slong i, slong j)
{
    /* TODO remove this slow debugging code */
    /*
    slong nrows = breadcrumb_mat_nrows(mat);
    slong ncols = breadcrumb_mat_ncols(mat);
    if (i < 0 || i >= nrows || j < 0 || j >= ncols)
    {
        flint_printf("tableau srcentry fail\n");
        flint_printf("i:%wd j:%wd nrows:%wd ncols:%wd\n", i, j, nrows, ncols);
        abort();
    }
    */
    return mat->data + i * mat->ncols + j;
}

static __inline__ breadcrumb_ptr
breadcrumb_mat_entry(breadcrumb_mat_t mat, slong i, slong j)
{
    /* TODO remove this slow debugging code */
    /*
    slong nrows = breadcrumb_mat_nrows(mat);
    slong ncols = breadcrumb_mat_ncols(mat);
    if (i < 0 || i >= nrows || j < 0 || j >= ncols)
    {
        flint_printf("tableau entry fail\n");
        flint_printf("i:%wd j:%wd nrows:%wd ncols:%wd\n", i, j, nrows, ncols);
        abort();
    }
    */
    return mat->data + i * mat->ncols + j;
}



#ifdef __cplusplus
}
#endif

#endif
