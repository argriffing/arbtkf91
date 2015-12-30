/*
 * Count the number of optimal solutions given a mask,
 * using fmpz integers.
 *
 * The memory usage is controlled as follows.
 * Find the row with the most 'contenders'.
 * Allocate an fmpz array long enough to hold two sparse rows of fmpzs.
 * Track two full rows of node pointers, some of which may be NULL.
 */

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"

#include "breadcrumbs.h"
#include "count_solutions.h"

static slong _breadcrumb_mat_max_row_nnz(const breadcrumb_mat_t mat);

slong
_breadcrumb_mat_max_row_nnz(const breadcrumb_mat_t mat)
{
    slong i, j, best, curr;
    best = 0;
    curr = 0;
    for (i = 0; i < breadcrumb_mat_nrows(mat); i++)
    {
        curr = 0;
        for (j = 0; j < breadcrumb_mat_ncols(mat); j++)
        {
            if (*breadcrumb_mat_srcentry(mat, i, j) & CRUMB_CONTENDER)
            {
                curr++;
            }
        }
        if (curr > best)
        {
            best = curr;
        }
    }
    return best;
}


void
count_solutions(fmpz_t res, const breadcrumb_mat_t mat)
{
    slong nrows, ncols, max_row_nnz;
    slong i, j, k;
    fmpz * counts;
    fmpz ** refs;
    breadcrumb_t d;
    fmpz *top, *diag, *left;
    fmpz * p;

    nrows = breadcrumb_mat_nrows(mat);
    ncols = breadcrumb_mat_nrows(mat);
    max_row_nnz = _breadcrumb_mat_max_row_nnz(mat);

    counts = _fmpz_vec_init(2 * max_row_nnz);
    refs = malloc(2 * ncols * sizeof(fmpz *));

    for (i = 0; i < nrows; i++)
    {
        k = 0;
        for (j = 0; j < ncols; j++)
        {
            d = *breadcrumb_mat_srcentry(mat, i, j);

            p = NULL;
            if (d & CRUMB_CONTENDER)
            {
                p = counts + (i%2)*max_row_nnz + k;
                fmpz_zero(p);
                k++;
            }
            refs[(i%2)*ncols + j] = p;

            if (p)
            {
                if (i == 0 && j == 0)
                {
                    fmpz_one(p);
                }

                top = diag = left = NULL;
                if (i && (d & CRUMB_TOP))
                {
                    top = refs[((i-1) % 2)*ncols + j];
                    fmpz_add(p, p, top);
                }
                if (i && j && (d & CRUMB_DIAG))
                {
                    diag = refs[((i-1) % 2)*ncols + (j-1)];
                    fmpz_add(p, p, diag);
                }
                if (j && (d & CRUMB_LEFT))
                {
                    left = refs[(i % 2)*ncols + j-1];
                    fmpz_add(p, p, left);
                }
            }
        }
    }

    i = nrows - 1;
    j = ncols - 1;
    p = refs[(i%2)*ncols + j];
    if (!p)
    {
        flint_printf("the bottom right entry of the tableau is empty\n");
        abort();
    }

    fmpz_set(res, p);

    _fmpz_vec_clear(counts, 2 * max_row_nnz);
    free(refs);
}
