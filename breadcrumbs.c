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


void
breadcrumb_mat_get_mask(breadcrumb_mat_t mask,
        const breadcrumb_mat_t mat, breadcrumb_t pattern)
{
    /* this is a generalized traceback, more like the backward phase
     * of a forward-backward algorithm */

    slong i, j, nr, nc;
    breadcrumb_t crumb;

    nr = breadcrumb_mat_nrows(mat);
    nc = breadcrumb_mat_ncols(mat);

    /* initialize the mask so that only the lower right entry is "on" */
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            *breadcrumb_mat_srcentry(mask, i, j) &= ~pattern;
        }
    }
    *breadcrumb_mat_entry(mask, nr-1, nc-1) |= pattern;

    for (i = nr-1; i >= 0; i--)
    {
        for (j = nc-1; j >= 0; j--)
        {
            if (*breadcrumb_mat_srcentry(mask, i, j) & pattern)
            {
                crumb = *breadcrumb_mat_srcentry(mat, i, j);

                /* printf for debug */
                flint_printf("%wd %wd %d\n", i, j, crumb);

                if (crumb & CRUMB_TOP)
                {
                    *breadcrumb_mat_entry(mask, i-1, j) |= pattern;
                }
                if (crumb & CRUMB_DIAG)
                {
                    *breadcrumb_mat_entry(mask, i-1, j-1) |= pattern;
                }
                if (crumb & CRUMB_LEFT)
                {
                    *breadcrumb_mat_entry(mask, i, j-1) |= pattern;
                }
            }
        }
    }
}


void
breadcrumb_mat_fprint(FILE *stream, const breadcrumb_mat_t mat)
{
    slong i, j;
    breadcrumb_t crumb;
    for (i = 0; i < breadcrumb_mat_nrows(mat); i++)
    {
        for (j = 0; j < breadcrumb_mat_ncols(mat); j++)
        {
            crumb = *breadcrumb_mat_srcentry(mat, i, j);
            flint_fprintf(stream, "%wd %wd %d\n", i, j, crumb);
        }
    }
}

slong
breadcrumb_mat_nnz(const breadcrumb_mat_t mat, breadcrumb_t pattern)
{
    slong i, j, nnz;
    nnz = 0;
    for (i = 0; i < breadcrumb_mat_nrows(mat); i++)
    {
        for (j = 0; j < breadcrumb_mat_ncols(mat); j++)
        {
            if (*breadcrumb_mat_srcentry(mat, i, j) & pattern)
            {
                nnz++;
            }
        }
    }
    return nnz;
}
