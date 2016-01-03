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
breadcrumb_mat_get_alignment(char *sa, char *sb, slong *plen,
        breadcrumb_mat_t mat, const slong *A, const slong *B)
{
    /*
     * Do the traceback.
     * The character arrays sa and sb are assumed to have been
     * already allocated.
     */
    slong i, j;
    char ACGT[4] = "ACGT";
    slong n = mat->nrows + mat->ncols;
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
    *plen = len;
}


void
breadcrumb_mat_get_mask(breadcrumb_mat_t mask, const breadcrumb_mat_t mat)
{
    /* this is a generalized traceback, more like the backward phase
     * of a forward-backward algorithm */

    slong i, j, nr, nc;
    breadcrumb_t crumb;

    nr = breadcrumb_mat_nrows(mat);
    nc = breadcrumb_mat_ncols(mat);

    if (breadcrumb_mat_nrows(mask) != nr ||
        breadcrumb_mat_ncols(mask) != nc)
    {
        flint_printf("traceback table dimensions mismatch\n");
        abort();
    }

    /* initialize the mask so that only the lower right entry is "on" */
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            *breadcrumb_mat_srcentry(mask, i, j) &= ~CRUMB_WANT2;
            *breadcrumb_mat_srcentry(mask, i, j) &= ~CRUMB_WANT3;
            *breadcrumb_mat_srcentry(mask, i, j) &= ~CRUMB_CONTENDER;
        }
    }
    *breadcrumb_mat_entry(mask, nr-1, nc-1) |= CRUMB_CONTENDER;

    for (i = nr-1; i >= 0; i--)
    {
        for (j = nc-1; j >= 0; j--)
        {
            crumb = *breadcrumb_mat_srcentry(mat, i, j);

            if (*breadcrumb_mat_srcentry(mask, i, j) & CRUMB_CONTENDER)
            {
                if (crumb & CRUMB_TOP)
                {
                    *breadcrumb_mat_entry(mask, i-1, j) |= CRUMB_CONTENDER;
                }
                if (crumb & CRUMB_DIAG)
                {
                    *breadcrumb_mat_entry(mask, i-1, j-1) |= CRUMB_CONTENDER;
                }
                if (crumb & CRUMB_LEFT)
                {
                    *breadcrumb_mat_entry(mask, i, j-1) |= CRUMB_CONTENDER;
                    *breadcrumb_mat_entry(mask, i, j-1) |= CRUMB_WANT2;
                }
            }

            if (*breadcrumb_mat_srcentry(mask, i, j) & CRUMB_WANT3)
            {
                if (crumb & CRUMB_TOP)
                {
                    *breadcrumb_mat_entry(mask, i-1, j) |= CRUMB_WANT3;
                }
                if (crumb & CRUMB_DIAG)
                {
                    *breadcrumb_mat_entry(mask, i-1, j-1) |= CRUMB_WANT3;
                }
                if (crumb & CRUMB_LEFT)
                {
                    *breadcrumb_mat_entry(mask, i, j-1) |= CRUMB_WANT2;
                }
            }

            if (*breadcrumb_mat_srcentry(mask, i, j) & CRUMB_WANT2)
            {
                if (crumb & CRUMB_DIAG2)
                {
                    *breadcrumb_mat_entry(mask, i-1, j-1) |= CRUMB_WANT3;
                }
                if (crumb & CRUMB_LEFT2)
                {
                    *breadcrumb_mat_entry(mask, i, j-1) |= CRUMB_WANT2;
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

void
breadcrumb_mat_set(breadcrumb_mat_t mat, const breadcrumb_mat_t src)
{
    slong nrows, ncols;

    nrows = breadcrumb_mat_nrows(src);
    ncols = breadcrumb_mat_ncols(src);

    if (breadcrumb_mat_nrows(mat) != nrows ||
        breadcrumb_mat_ncols(mat) != ncols)
    {
        flint_printf("traceback table dimensions mismatch\n");
        abort();
    }

    memcpy(mat->data, src->data, nrows * ncols * sizeof(breadcrumb_t));
}
