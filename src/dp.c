#include <stdlib.h>
#include <string.h>

#include "flint/flint.h"

#include "dp.h"


void
dp_mat_init(dp_mat_t mat, slong nrows, slong ncols)
{
    /*
     * Initially we are interested in max2 and max3 for every cell.
     * Each of {m1, m2} is a candidate for max2,
     * and each of {m0, m1, m2} is a candidate for max3.
     * All tableau cells are possible trace candidates.
     */
    slong n = nrows * ncols;
    mat->data = malloc(n * sizeof(dp_t));
    mat->nrows = nrows;
    mat->ncols = ncols;
    int i;
    for (i = 0; i < n; i++)
    {
        mat->data[i] = 0xFF;
    }
}

void
dp_mat_clear(dp_mat_t mat)
{
    free(mat->data);
}

void
dp_mat_get_alignment(char *sa, char *sb, slong *plen,
        dp_mat_t mat, const slong *A, const slong *B)
{
    /*
     * Do the traceback. The character arrays sa and sb
     * are assumed to have been already allocated.
     */
    slong i, j;
    char ACGT[4] = "ACGT";
    slong len = 0;
    i = mat->nrows - 1;
    j = mat->ncols - 1;
    dp_t x;
    while (i > 0 || j > 0)
    {
        x = *dp_mat_entry(mat, i, j);
        if (x & DP_MAX3_M0)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = '-';
            i--;
        }
        else if (x & DP_MAX3_M1)
        {
            sa[len] = ACGT[A[i-1]];
            sb[len] = ACGT[B[j-1]];
            i--;
            j--;
        }
        else if (x & DP_MAX3_M2)
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
dp_mat_backward(dp_mat_t mat)
{
    /* This is a generalized traceback. */
    slong i, j, nr, nc;
    dp_t *px;
    dp_t y;

    nr = dp_mat_nrows(mat);
    nc = dp_mat_ncols(mat);

    /* Clear the flags that will be updated in this function. */
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            px = dp_mat_entry(mat, i, j);
            *px &= ~(DP_MAX3 | DP_TRACE | DP_MAX2);
        }
    }

    /* Set the flags for the bottom right cell. */
    *dp_mat_entry(mat, nr-1, nc-1) |= (DP_MAX3 | DP_TRACE);

    /* Update flags recursively. */
    for (i = nr-1; i >= 0; i--)
    {
        for (j = nc-1; j >= 0; j--)
        {
            px = dp_mat_entry(mat, i, j);

            /* Check the cell to the right, if available. */
            if (j < nc-1)
            {
                y = *dp_mat_entry(mat, i, j+1);
                if ((y & DP_TRACE) && (y & DP_MAX3_M2))
                {
                    *px |= DP_TRACE;
                }
                if (dp_m2_is_interesting(y))
                {
                    *px |= DP_MAX2;
                }
            }

            /* Check the cell to the lower right, if available. */
            if (i < nr-1 && j < nc-1)
            {
                y = *dp_mat_entry(mat, i+1, j+1);
                if ((y & DP_TRACE) && (y & DP_MAX3_M1))
                {
                    *px |= DP_TRACE;
                }
                if (dp_m1_is_interesting(y))
                {
                    *px |= DP_MAX3;
                }
            }

            /* Check the cell below, if available. */
            if (i < nr-1)
            {
                y = *dp_mat_entry(mat, i+1, j);
                if ((y & DP_TRACE) && (y & DP_MAX3_M0))
                {
                    *px |= DP_TRACE;
                }
                if (dp_m0_is_interesting(y))
                {
                    *px |= DP_MAX3;
                }
            }
        }
    }
}


void
dp_mat_set(dp_mat_t mat, const dp_mat_t src)
{
    slong nrows, ncols;

    nrows = dp_mat_nrows(src);
    ncols = dp_mat_ncols(src);

    if (dp_mat_nrows(mat) != nrows ||
        dp_mat_ncols(mat) != ncols)
    {
        flint_printf("traceback table dimensions mismatch\n");
        abort();
    }

    memcpy(mat->data, src->data, nrows * ncols * sizeof(dp_t));
}

void
dp_mat_check_alignment(
        int *p_is_optimal, int *p_is_canonical,
        dp_mat_t mat, const slong *A, const slong *B, slong len)
{
    slong i, j, k;
    i = mat->nrows - 1;
    j = mat->ncols - 1;
    dp_t full, canonical, observed;

    *p_is_optimal = 1;
    *p_is_canonical = 1;

    k = len-1;
    while (i > 0 || j > 0)
    {
        observed = 0;
        if (A[k] > -1 && B[k] == -1)
        {
            observed = DP_MAX3_M0;
        }
        else if (A[k] > -1 && B[k] > -1)
        {
            observed = DP_MAX3_M1;
        }
        else if (A[k] == -1 && B[k] > -1)
        {
            observed = DP_MAX3_M2;
        }
        else
        {
            flint_printf("unexpected alignment column\n");
            abort();
        }

        full = *dp_mat_entry(mat, i, j);
        if (full & DP_MAX3_M0)
        {
            canonical = DP_MAX3_M0;
        }
        else if (full & DP_MAX3_M1)
        {
            canonical = DP_MAX3_M1;
        }
        else if (full & DP_MAX3_M2)
        {
            canonical = DP_MAX3_M2;
        }
        else
        {
            flint_printf("lost the thread ");
            flint_printf("in the dynamic programing traceback\n");
            abort();
        }

        if (observed != canonical)
        {
            *p_is_canonical = 0;
        }
        if (!(observed & full))
        {
            *p_is_optimal = 0;
            return;
        }

        if (observed == DP_MAX3_M0)
        {
            i--;
        }
        else if (observed == DP_MAX3_M1)
        {
            i--;
            j--;
        }
        else if (observed == DP_MAX3_M2)
        {
            j--;
        }
        k--;
    }
}
