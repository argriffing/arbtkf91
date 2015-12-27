/*
 * Alignment disambiguation based on a sparse (CSR) representation
 * of possible cells of a dynamic programming table, after having
 * used upper and lower probability bounds to determine that the
 * maximum probability traceback must pass through these cells.
 */

#include "flint/flint.h"

#include "breadcrumbs.h"


typedef struct bound_node_struct
{
    fmpz * pmax2;
    fmpz * pmax3;
    bound_node_struct * top;
    bound_node_struct * diag;
    bound_node_struct * left;

    /* used to indicate whether the optimal traceback may
     * possibly contain this node */
    int mark; 

} bound_node_struct;
typedef bound_node_struct * bound_node_ptr;
typedef bound_node_struct bound_node_t[1];

void
bound_node_init(bound_node_t x)
{
    x->pmax2 = NULL;
    x->pmax3 = NULL;
    x->top = NULL;
    x->diag = NULL;
    x->left = NULL;
    x->mark = 0;
}


typedef struct
{
    slong * indptr; /* length r + 1 */
    slong * indices; /* length nnz */
    bound_node_struct * data; /* length nnz */
    fmpz * brick_of_fmpz; /* length 2*2*rank*max_row_nnz */
    slong r; /* number of rows */
    slong c; /* number of columns */
    slong nnz; /* number of nonzero entries */
    slong max_row_nnz; /* max number of nonzero entries in a single row */
    slong rank; /* length of an fmpz vector that tracks partial aln state */
} bound_mat_struct;
typedef bound_mat_struct * bound_mat_ptr;
typedef bound_mat_struct bound_mat_t[1];

slong
_breadcrumb_mat_max_row_nnz(const breadcrumb_mat_t mask, breadcrumb_t pattern)
{
    slong i, j, best, curr;
    best = 0;
    curr = 0;
    for (i = 0; i < breadcrumb_mat_nrows(mat); i++)
    {
        curr = 0;
        for (j = 0; j < breadcrumb_mat_ncols(mat); j++)
        {
            if (*breadcrumb_mat_srcentry(mat, i, j) & pattern)
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


fmpz_ptr
_bound_mat_max2ref(bound_mat_t mat, i, k)
{
    /* this function is intended to be called only during initialization */
    /* i : row index of the entry */
    /* k : the entry follows k nonzero entries in row i */
    fmpz_ptr base;
    base = mat->brick_of_fmpz + (i % 2) * mat->max_row_nnz * 2 * rank;
    return base + k * 2 * rank
}

fmpz_ptr
_bound_mat_max3ref(bound_mat_t mat, i, k)
{
    /* this function is intended to be called only during initialization */
    /* i : row index of the entry */
    /* k : the entry follows k nonzero entries in row i */
    fmpz_ptr base;
    base = mat->brick_of_fmpz + (i % 2) * mat->max_row_nnz * 2 * rank;
    return base + k * 2 * rank + rank;
}

void
bound_mat_init(bound_mat_t mat, breadcrumb_mat_t mask, breadcrumb_t pattern,
        slong rank)
{
    /* pattern : a bit pattern that indicates whether a cell is potentially
     *          part of the maximum likelihood traceback
     * rank : length of fmpz vector that tracks partial alignment state
     */
    slong nr, nc, nnz, max_row_nnz;
    slong i, j, k;
    slong idx;
    breadcrumb_t crumb;
    bound_node_ptr p;
    bound_node_ptr tmp;

    nr = breadcrumb_mat_nrows(mask);
    nc = breadcrumb_mat_ncols(mask);
    nnz = breadcrumb_mat_nnz(mask, pattern);
    max_row_nnz = _breadcrumb_mat_max_row_nnz(mask, pattern);

    mat->r = nr;
    mat->c = nc;
    mat->nnz = nnz;
    mat->rank = rank;
    mat->max_row_nnz = max_row_nnz;

    /* allocate arrays */
    mat->indptr = flint_malloc((nr + 1) * sizeof(slong));
    mat->indices = flint_malloc(nnz * sizeof(slong));
    mat->data = flint_malloc(nnz * sizeof(bound_node_struct));
    mat->brick_of_fmpz = _fmpz_vec_init(2 * max_row_nnz * 2 * rank);

    /* temporary array, to help link together the nodes */
    tmp = flint_malloc(nc * 2 * sizeof(bound_node_ptr));

    /* init the data */
    for (i = 0; i < nnz; i++)
    {
        bound_node_init(mat->data + i);
    }

    /* initialize the CSR matrix */
    idx = 0;
    mat->indptr[0] = 0;
    for (i = 0; i < nr; i++)
    {
        k = 0;
        for (j = 0; j < nc; j++)
        {
            crumb = *breadcrumb_mat_srcentry(mask, i, j);
            if (crumb & pattern)
            {
                p = mat->data + idx;
                p->pmat2 = _bound_mat_max2ref(mat, i, k)
                p->pmat3 = _bound_mat_max3ref(mat, i, k)
                tmp[(i % 2) * nc + j] = p;
                mat->indices[idx] = j;
                if (crumb & CRUMB_TOP)
                {
                    p->top = tmp + ((i-1) % 2) * nc + j;
                }
                if (crumb & CRUMB_DIAG)
                {
                    p->diag = tmp + ((i-1) % 2) * nc + (j-1);
                }
                if (crumb & CRUMB_LEFT)
                {
                    p->left = tmp + (i % 2) * nc + (j-1);
                }
                idx++;
                k++;
            }
        }
        mat->indptr[i+1] = idx;
    }

    flint_free(tmp);
}

void
bound_mat_init(bound_mat_t mat)
{
}

void
bound_mat_clear(bound_mat_t mat)
{
    /* free arrays */
    flint_free(mat->indptr);
    flint_free(mat->indices);
    flint_free(mat->data);
    _fmpz_vec_clear(mat->brick_of_fmpz);
}
