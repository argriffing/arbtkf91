/*
 * Alignment disambiguation based on a sparse (CSR) representation
 * of possible cells of a dynamic programming table, after having
 * used upper and lower probability bounds to determine that the
 * maximum probability traceback must pass through these cells.
 */

#include "flint/flint.h"

#include "breadcrumbs.h"
#include "tkf91_generator_vecs.h"
#include "bound_mat.h"

slong _breadcrumb_mat_nnz(const breadcrumb_mat_t mat);
slong _breadcrumb_mat_max_row_nnz(const breadcrumb_mat_t mat);

typedef struct bound_node_struct_tag
{
    fmpz * pmax2;
    fmpz * pmax3;
    struct bound_node_struct_tag * top;
    struct bound_node_struct_tag * diag;
    struct bound_node_struct_tag * left;
    breadcrumb_t crumb;
} bound_node_struct;
typedef bound_node_struct * bound_node_ptr;
typedef bound_node_struct bound_node_t[1];

void bound_node_init(bound_node_t x);

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

void bound_mat_init(bound_mat_t mat, breadcrumb_mat_t mask, slong rank);
void bound_mat_clear(bound_mat_t mat);
fmpz * _bound_mat_max3ref(bound_mat_t mat, slong i, slong k);
fmpz * _bound_mat_max2ref(bound_mat_t mat, slong i, slong k);
slong bound_mat_ncols(const bound_mat_t mat);
slong bound_mat_nrows(const bound_mat_t mat);

void _check_equal(fmpz * u, fmpz * v, slong r);
void _tkf91_dp_verify_symbolically(
        tkf91_generator_vecs_t h,
        bound_mat_t b,
        const slong *A,
        const slong *B);

slong
bound_mat_nrows(const bound_mat_t mat)
{
    return mat->r;
}

slong
bound_mat_ncols(const bound_mat_t mat)
{
    return mat->c;
}

void
bound_node_init(bound_node_t x)
{
    x->pmax2 = NULL;
    x->pmax3 = NULL;
    x->top = NULL;
    x->diag = NULL;
    x->left = NULL;
    x->crumb = 0;
}

slong
_breadcrumb_mat_nnz(const breadcrumb_mat_t mat)
{
    breadcrumb_t pattern;
    slong i, j, nnz;
    pattern = CRUMB_WANT2 | CRUMB_WANT3 | CRUMB_CONTENDER;
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


slong
_breadcrumb_mat_max_row_nnz(const breadcrumb_mat_t mat)
{
    breadcrumb_t pattern;
    slong i, j, best, curr;
    pattern = CRUMB_WANT2 | CRUMB_WANT3 | CRUMB_CONTENDER;
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

fmpz *
_bound_mat_max2ref(bound_mat_t mat, slong i, slong k)
{
    /* this function is intended to be called only during initialization */
    /* i : row index of the entry */
    /* k : the entry follows k nonzero entries in row i */
    slong rank;
    fmpz * base;
    rank = mat->rank;
    base = mat->brick_of_fmpz + (i % 2) * mat->max_row_nnz * 2 * rank;
    return base + k * 2 * rank;
}

fmpz *
_bound_mat_max3ref(bound_mat_t mat, slong i, slong k)
{
    /* this function is intended to be called only during initialization */
    /* i : row index of the entry */
    /* k : the entry follows k nonzero entries in row i */
    slong rank;
    fmpz * base;
    rank = mat->rank;
    base = mat->brick_of_fmpz + (i % 2) * mat->max_row_nnz * 2 * rank;
    return base + k * 2 * rank + rank;
}

void
bound_mat_init(bound_mat_t mat, breadcrumb_mat_t mask, slong rank)
{
    slong nr, nc, nnz, max_row_nnz;
    slong i, j, k;
    slong idx;
    breadcrumb_t crumb;
    bound_node_ptr p;
    bound_node_ptr * tmp;

    nr = breadcrumb_mat_nrows(mask);
    nc = breadcrumb_mat_ncols(mask);
    nnz = _breadcrumb_mat_nnz(mask);
    max_row_nnz = _breadcrumb_mat_max_row_nnz(mask);

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
            tmp[(i % 2) * nc + j] = NULL;
        }
        for (j = 0; j < nc; j++)
        {
            crumb = *breadcrumb_mat_srcentry(mask, i, j);
            if (crumb & (CRUMB_CONTENDER | CRUMB_WANT3 | CRUMB_WANT2))
            {
                mat->indices[idx] = j;
                p = mat->data + idx;
                tmp[(i % 2) * nc + j] = p;
                p->pmax2 = _bound_mat_max2ref(mat, i, k);
                p->pmax3 = _bound_mat_max3ref(mat, i, k);
                if (0 < i)
                {
                    p->top = tmp[((i-1) % 2) * nc + j];
                }
                if (0 < i && 0 < j)
                {
                    p->diag = tmp[((i-1) % 2) * nc + (j-1)];
                }
                if (0 < j)
                {
                    p->left = tmp[(i % 2) * nc + (j-1)];
                }
                p->crumb = crumb;
                idx++;
                k++;
            }
        }
        mat->indptr[i+1] = idx;
    }

    flint_free(tmp);
}


void
bound_mat_clear(bound_mat_t mat)
{
    /* free arrays */
    flint_free(mat->indptr);
    flint_free(mat->indices);
    flint_free(mat->data);
    _fmpz_vec_clear(mat->brick_of_fmpz,
            2 * mat->max_row_nnz * 2 * mat->rank);
}


void
tkf91_dp_verify_symbolically(
        fmpz_mat_t mat,
        const tkf91_generator_indices_t g,
        breadcrumb_mat_t mask,
        const slong *A,
        const slong *B)
{
    /* Inputs:
     *   mat : the generator matrix -- mat_ij where i is a generator index
     *         and j is an expression index.
     *   g : a struct with tkf91 generator indices
     *   mask : a previously computed generalized traceback,
     *          with marks indicating whether cells are possibly in the
     *          max likelihood traceback, and with directional links
     *          indicating which direction(s) are best, backwards,
     *          from each cell.
     */
    fmpz_mat_t H, V;
    slong rank;
    tkf91_generator_vecs_t h;
    bound_mat_t b;

    /* Compute a Hermite decomposition of the generator matrix. */
    /* U*mat = H ; U^-1 = V ; rank = rank(H) */
    fmpz_mat_init(H, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
    fmpz_mat_init(V, fmpz_mat_nrows(mat), fmpz_mat_nrows(mat));
    _fmpz_mat_hnf_inverse_transform(H, V, &rank, mat);

    /*
    flint_printf("matrix rank ");
    flint_printf("determined from the Hermite normal form: %wd\n", rank);
    */

    tkf91_generator_vecs_init(h, g, V, rank);
    bound_mat_init(b, mask, rank);

    _tkf91_dp_verify_symbolically(h, b, A, B);

    fmpz_mat_clear(H);
    fmpz_mat_clear(V);
    tkf91_generator_vecs_clear(h);
    bound_mat_clear(b);
}


void
_check_equal(fmpz * u, fmpz * v, slong r)
{
    /*
    flint_printf("verifying... ");
    */
    if (!_fmpz_vec_equal(u, v, r))
    {
        flint_printf("verification failed\n");
        abort();
    }
    /*
    flint_printf("OK\n");
    */
}


void
_tkf91_dp_verify_symbolically(
        tkf91_generator_vecs_t h,
        bound_mat_t b,
        const slong *A,
        const slong *B)
{
    /*
     * Use only a forward pass, not a backward pass.
     * This should be enough to check whether
     * the directions indicated as plausible by the nodes of the bound_mat
     * actually correspond to identical integer state vectors.
     */
    slong r;
    slong nr;

    bound_node_ptr p0, p1, p2;

    fmpz * m0;
    fmpz * m1;
    fmpz * m2;

    fmpz * m0_buf;
    fmpz * m1_buf;
    fmpz * m2_buf;

    r = tkf91_generator_vecs_rank(h);
    nr = bound_mat_nrows(b);

    m0_buf = _fmpz_vec_init(r);
    m1_buf = _fmpz_vec_init(r);
    m2_buf = _fmpz_vec_init(r);

    /* iterate over rows of the CSR dynamic programming matrix */
    slong i, j, u;
    slong nta, ntb;
    slong start, stop;
    bound_node_ptr node;
    breadcrumb_t crumb, want2, want3;
    for (i = 0; i < nr; i++)
    {
        start = b->indptr[i];
        stop = b->indptr[i+1];
        for (u = start; u < stop; u++)
        {
            j = b->indices[u];
            node = b->data+u;

            crumb = node->crumb;
            want2 = crumb & CRUMB_WANT2;
            want3 = crumb & (CRUMB_CONTENDER | CRUMB_WANT3);

            /*
             * At this point, the vector pointers point to an appropriate
             * location in a cyclic buffer.
             * For each of the two vectors, if the vector is requested
             * then initialize its entries to zero, otherwise set
             * the vector pointer to NULL.
             */
            if (want2)
            {
                _fmpz_vec_zero(node->pmax2, r);
            }
            else
            {
                node->pmax2 = NULL;
            }
            if (want3)
            {
                _fmpz_vec_zero(node->pmax3, r);
            }
            else
            {
                node->pmax3 = NULL;
            }

            /*
             * Set m0, m1, m2 to either the appropriate utility fmpz vec
             * or to NULL, depending upon whether the vector is requested.
             */
            m0 = m1 = m2 = NULL;
            if (want3 && (node->crumb & CRUMB_TOP))
            {
                m0 = m0_buf;
            }
            if ((want3 && (node->crumb & CRUMB_DIAG)) ||
                (want2 && (node->crumb & CRUMB_DIAG2)))
            {
                m1 = m1_buf;
            }
            if ((want3 && (node->crumb & CRUMB_LEFT)) ||
                (want2 && (node->crumb & CRUMB_LEFT2)))
            {
                m2 = m2_buf;
            }

            if (i < 1 || j < 1)
            {
                if (i == 0 && j == 0)
                {
                    if (m1) _fmpz_vec_set(m1, h->m1_00, r);
                }
                else if (i == 1 && j == 0)
                {
                    if (m0) _fmpz_vec_set(m0, h->m0_10, r);
                }
                else if (i == 0 && j == 1)
                {
                    if (m2) _fmpz_vec_set(m2, h->m2_01, r);
                }
                else
                {
                    if (i == 0)
                    {
                        if (m2)
                        {
                            ntb = B[j - 1];
                            p2 = node->left;
                            _fmpz_vec_add(m2, p2->pmax2, h->m2_0j_incr[ntb], r);
                        }
                    }
                    else if (j == 0)
                    {
                        if (m0)
                        {
                            nta = A[i - 1];
                            p0 = node->top;
                            _fmpz_vec_add(m0, p0->pmax3, h->m0_i0_incr[nta], r);
                        }
                    }
                }
            }
            else
            {
                nta = A[i - 1];
                ntb = B[j - 1];

                if (m0)
                {
                    p0 = node->top;
                    _fmpz_vec_add(m0, p0->pmax3, h->c0_incr[nta], r);
                }

                if (m1)
                {
                    p1 = node->diag;
                    _fmpz_vec_add(m1, p1->pmax3, h->c1_incr[nta*4+ntb], r);
                }

                if (m2)
                {
                    p2 = node->left;
                    _fmpz_vec_add(m2, p2->pmax2, h->c2_incr[ntb], r);
                }
            }

            /*
             * This is the verification step.
             * In another module, numerical ties among subsets of
             * {m0, m1, m2} may have been detected, and these ties
             * will be checked here 'symbolically' using integer vectors.
             * If elements of subsets that are numerically guessed
             * to be identical have identical corresponding integer vectors,
             * then the verification succeeds.
             * Otherwise the verification fails.
             * There are multiple possible explanations for verification
             * failure.
             * 1) a bug
             * 2) a true numerical 'coincidence' slipped through, because
             *    the factor refinement is only clever enough to work
             *    with rational numbers but not sparse polynomials
             * 3) the numerical bounds were not computed with
             *    high enough precision
             * Ideally the explanation is (3), but the other two explanations
             * are also possible.
             */
            if (want2)
            {
                if ((crumb & CRUMB_DIAG2) && (crumb & CRUMB_LEFT2))
                {
                    _check_equal(m1, m2, r);
                }
            }
            if (want3)
            {
                if ((crumb & CRUMB_TOP) && (crumb & CRUMB_DIAG))
                {
                    _check_equal(m0, m1, r);
                }
                if ((crumb & CRUMB_DIAG) && (crumb & CRUMB_LEFT))
                {
                    _check_equal(m1, m2, r);
                }
                if ((crumb & CRUMB_LEFT) && (crumb & CRUMB_TOP))
                {
                    _check_equal(m2, m0, r);
                }
            }

            /*
             * Update the pmax2 and pmax3 vectors.
             * It is not required to take a maximum because we just
             * finished verifying equality of maxima.
             */
            if (want2)
            {
                if (m1)
                {
                    _fmpz_vec_set(node->pmax2, m1, r);
                }
                else if (m2)
                {
                    _fmpz_vec_set(node->pmax2, m2, r);
                }
            }
            if (want3)
            {
                if (m0)
                {
                    _fmpz_vec_set(node->pmax3, m0, r);
                }
                else if (m1)
                {
                    _fmpz_vec_set(node->pmax3, m1, r);
                }
                else if (m2)
                {
                    _fmpz_vec_set(node->pmax3, m2, r);
                }
            }
        }
    }

    _fmpz_vec_clear(m0_buf, r);
    _fmpz_vec_clear(m1_buf, r);
    _fmpz_vec_clear(m2_buf, r);
}
