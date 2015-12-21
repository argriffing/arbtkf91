#include <stdlib.h>

#include "flint/flint.h"

#include "breadcrumbs.h"


void
breadcrumb_mat_init(breadcrumb_mat_t mat, slong nrows, slong ncols)
{
    mat->data = calloc(nrows * ncols, sizeof(breadcrumb_t));
    mat->nrows = nrows;
    mat->ncols = ncols;
}

void
breadcrumb_mat_clear(breadcrumb_mat_t mat)
{
    free(mat->data);
}

breadcrumb_ptr
breadcrumb_mat_entry(breadcrumb_mat_t mat, slong i, slong j)
{
    /*
    if (i < 0 || i >= mat->nrows || j < 0 || j >= mat->ncols)
    {
        flint_printf("breadcrumb matrix indexing error\n");
        flint_printf("i=%wd j=%wd\n", i, j);
        flint_printf("nrows=%wd ncols=%wd\n", mat->nrows, mat->ncols);
        abort();
    }
    */
    return mat->data + i * mat->ncols + j;
}

void
breadcrumb_mat_get_alignment(char **psa, char **psb,
        breadcrumb_mat_t mat, const slong *A, const slong *B)
{
    /* do the traceback */
    slong i, j;
    char ACGT[4] = "ACGT";
    slong n = mat->nrows + mat->ncols;
    char *sa = calloc(n, sizeof(char));
    char *sb = calloc(n, sizeof(char));
    slong len = 0;
    i = mat->nrows - 1;
    j = mat->ncols - 1;
    breadcrumb_t crumb;
    while (i > 0 || j > 0)
    {
        crumb = *breadcrumb_mat_entry(mat, i, j);
        if (crumb & CRUMB_TOP)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = '-';
            i--;
        }
        else if (crumb & CRUMB_DIAG)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = ACGT[B[j-1]];
            i--;
            j--;
        }
        else if (crumb & CRUMB_LEFT)
        {
            sa[len] = '-';
            sb[len] = ACGT[B[j-1]];
            j--;
        }
        else
        {
            flint_printf("lost the thread ");
            flint_printf("in the dynamic programing traceback\n");
            abort();
        }
        len++;
    }
    char tmp;
    for (i = 0; i < len/2; i++)
    {
        j = len - 1 - i;
        tmp = sa[i]; sa[i] = sa[j]; sa[j] = tmp;
        tmp = sb[i]; sb[i] = sb[j]; sb[j] = tmp;
    }
    *psa = sa;
    *psb = sb;
}
