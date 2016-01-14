/*
 * Arbitrary precision tkf91 dynamic programming.
 *
 * After some experimentation, this method of using arb_t values has been
 * shown to be slower than explicitly tracking upper and lower mag_t bounds.
 */

#include <time.h>

#include "arb_mat.h"

#include "tkf91_dp.h"
#include "tkf91_dp_r.h"
#include "dp.h"
#include "wavefront_hermite.h"
#include "tkf91_generator_vecs.h"
#include "unused.h"




void
_hwave_element_add_vec(
        hwave_element_t e, const fmpz * vec1, const fmpz * vec2, arb_srcptr v,
        slong r, slong prec);

void
_hwave_element_set_vec(
        hwave_element_t e, const fmpz * vec, arb_srcptr v,
        slong r, slong prec);

void
_hwave_element_add_vec(
        hwave_element_t e, const fmpz * vec1, const fmpz * vec2, arb_srcptr v,
        slong r, slong prec)
{
    _fmpz_vec_add(e->vec, vec1, vec2, r);
    _arb_vec_dot_fmpz_vec(e->value, v, e->vec, r, prec);
    e->status = HWAVE_STATUS_UNAMBIGUOUS;
}

void
_hwave_element_set_vec(
        hwave_element_t e, const fmpz * vec, arb_srcptr v,
        slong r, slong prec)
{
    /* Input:
     *  e : one of the three elements of an hwave cell
     *  vec : integer vector of r coefficients
     *  v : arb vector of r real balls
     *  r : the rank of the system
     *  prec : precision for arb
     */

    /* copy the integer coefficients */
    _fmpz_vec_set(e->vec, vec, r);

    /* update the dot product */
    _arb_vec_dot_fmpz_vec(e->value, v, e->vec, r, prec);

    /*
     * Mark that for this element the vector of integer coefficients
     * is known unambiguously.
     */
    e->status = HWAVE_STATUS_UNAMBIGUOUS;
}

hwave_element_ptr
_get_max3_checked(
        hwave_element_ptr e0,
        hwave_element_ptr e1,
        hwave_element_ptr e2,
        slong len);

hwave_element_ptr
_get_max3_checked(
        hwave_element_ptr e0,
        hwave_element_ptr e1,
        hwave_element_ptr e2,
        slong len)
{
    /*
     * FIXME
     * For now we abort if an ambiguity is detected numerically
     * but not symbolically.
     */
    /* find the element with the greatest midpoint */
    slong i;
    hwave_element_ptr x[] = {e0, e1, e2};

    /* for now require that all elements have unambiguous status */
    for (i = 0; i < 3; i++)
    {
        if (x[i]->status == HWAVE_STATUS_AMBIGUOUS)
        {
            flint_printf("ambiguous status -- ");
            flint_printf("handling this case is not yet implemented...\n");
            abort();
        }
    }

    /* determine the element with greatest midpoint */
    hwave_element_ptr best;
    best = e0;
    if (arf_cmp(arb_midref(best->value), arb_midref(e1->value)) < 0) {
        best = e1;
    }
    if (arf_cmp(arb_midref(best->value), arb_midref(e2->value)) < 0) {
        best = e2;
    }

    /*
     * If any element's value overlaps the value of the element
     * whose midpoint is greatest, then complain if the integer coefficient
     * vectors are not equal.
     */
    for (i = 0; i < 3; i++)
    {
        if (arb_overlaps(best->value, x[i]->value))
        {
            if (!_fmpz_vec_equal(best->vec, x[i]->vec, len))
            {
                flint_printf("an ambiguity was detected numericaly ");
                flint_printf("but not symbolically -- ");
                flint_printf("dealing with this situation ");
                flint_printf("is not yet implemented...\n");
                abort();
            }
        }
    }

    return best;
}

void tkf91_dynamic_programming_hermite(
        solution_t sol, const request_t req,
        tkf91_generator_vecs_t g,
        arb_ptr v, slong rank, slong prec,
        const slong *A, slong szA,
        const slong *B, slong szB);

void tkf91_dynamic_programming_hermite(
        solution_t sol, const request_t req,
        tkf91_generator_vecs_t g,
        arb_ptr v, slong rank, slong prec,
        const slong *A, slong szA,
        const slong *B, slong szB)
{
    UNUSED(req);

    /*
     * Define the matrix to be used for the traceback.
     * The number of rows is one greater than the length of the first sequence,
     * and the number of columns is one greater than the length of the second
     * sequence.
     */
    slong nrows = szA + 1;
    slong ncols = szB + 1;
    dp_mat_t crumb_mat;
    dp_mat_init(crumb_mat, nrows, ncols);

    /* define the wavefront matrix */
    hwave_mat_t wave;
    slong modulus = 3;
    /* slong modulus = nrows + ncols - 1; */
    hwave_mat_init(wave, nrows + ncols - 1, modulus, rank);

    /* iterate over anti-diagonal bands of the dynamic programming matrix */
    hwave_element_ptr best;
    hwave_cell_ptr cell, p0, p1, p2;
    slong i, j, k, l;
    slong istart, jstart, lstart;
    for (k = 0; k < nrows + ncols - 1; k++)
    {
        if (k < nrows)
        {
            istart = k;
            jstart = 0;
        }
        else
        {
            istart = nrows - 1;
            jstart = k - (nrows - 1);
        }
        lstart = nrows - 1 + jstart - istart;

        /* iterate over entries of the diagonal */
        i = istart;
        j = jstart;
        l = lstart;
        slong nta, ntb;
        while (0 <= i && i < nrows && 0 <= j && j < ncols)
        {
            /* check some invariants */
            if (k != i + j)
            {
                flint_printf("wavefront indexing problem ");
                flint_printf("i=%wd j=%wd k=%wd\n", i, j, k);
                abort();
            }
            if (l != nrows - 1 + j - i)
            {
                flint_printf("wavefront indexing problem ");
                flint_printf("i=%wd j=%wd l=%wd\n", i, j, l);
                abort();
            }

            cell = hwave_mat_entry(wave, k, l);
            if (i == 0 && j == 0)
            {
                hwave_element_set_undefined(cell->m+0);
                _hwave_element_set_vec(cell->m+1, g->m1_00, v, rank, prec);
                hwave_element_set_undefined(cell->m+2);
            }
            else if (i == 1 && j == 0)
            {
                _hwave_element_set_vec(cell->m+0, g->m0_10, v, rank, prec);
                hwave_element_set_undefined(cell->m+1);
                hwave_element_set_undefined(cell->m+2);
            }
            else if (i == 0 && j == 1)
            {
                hwave_element_set_undefined(cell->m+0);
                hwave_element_set_undefined(cell->m+1);
                _hwave_element_set_vec(cell->m+2, g->m2_01, v, rank, prec);
            }
            else
            {
                if (i == 0)
                {
                    ntb = B[j - 1];
                    p2 = hwave_mat_entry_left(wave, k, l);
                    hwave_element_set_undefined(cell->m+0);
                    hwave_element_set_undefined(cell->m+1);
                    if (p2->m[2].status != HWAVE_STATUS_UNAMBIGUOUS)
                    {
                        flint_printf("found a not unambiguous status ");
                        flint_printf("while filling the first row\n");
                        abort();
                    }
                    _hwave_element_add_vec(
                            cell->m+2,
                            p2->m[2].vec,
                            g->m2_0j_incr[ntb],
                            v, rank, prec);
                }
                else if (j == 0)
                {
                    nta = A[i - 1];
                    p0 = hwave_mat_entry_top(wave, k, l);
                    if (p0->m[0].status != HWAVE_STATUS_UNAMBIGUOUS)
                    {
                        flint_printf("found a not unambiguous status ");
                        flint_printf("while filling the first column\n");
                        abort();
                    }
                    _hwave_element_add_vec(
                            cell->m+0,
                            p0->m[0].vec,
                            g->m0_i0_incr[nta],
                            v, rank, prec);
                    hwave_element_set_undefined(cell->m+1);
                    hwave_element_set_undefined(cell->m+2);
                }
                else
                {
                    nta = A[i - 1];
                    ntb = B[j - 1];
                    p0 = hwave_mat_entry_top(wave, k, l);
                    p1 = hwave_mat_entry_diag(wave, k, l);
                    p2 = hwave_mat_entry_left(wave, k, l);

                    best = _get_max3_checked(p0->m+0, p0->m+1, p0->m+2, rank);
                    if (best->status == HWAVE_STATUS_UNDEFINED)
                    {
                        hwave_element_set_undefined(cell->m+0);
                    }
                    else
                    {
                        _hwave_element_add_vec(
                                cell->m+0,
                                best->vec,
                                g->c0_incr[nta],
                                v, rank, prec);
                    }

                    best = _get_max3_checked(p1->m+0, p1->m+1, p1->m+2, rank);
                    if (best->status == HWAVE_STATUS_UNDEFINED)
                    {
                        hwave_element_set_undefined(cell->m+1);
                    }
                    else
                    {
                        _hwave_element_add_vec(
                                cell->m+1,
                                best->vec,
                                g->c1_incr[nta*4 + ntb],
                                v, rank, prec);
                    }

                    /* use max3 to check only max2 by duplicating one -- */
                    /* passing p2->m+1 twice is not a bug */
                    /* TODO make a max2 convenience function for this... */
                    best = _get_max3_checked(p2->m+1, p2->m+1, p2->m+2, rank);
                    if (best->status == HWAVE_STATUS_UNDEFINED)
                    {
                        hwave_element_set_undefined(cell->m+2);
                    }
                    else
                    {
                        _hwave_element_add_vec(
                                cell->m+2,
                                best->vec,
                                g->c2_incr[ntb],
                                v, rank, prec);
                    }
                }
            }

            /* fill the table for traceback */
            /* FIXME use the breadcrumb ambiguity bit as appropriate */
            /* FIXME This entire module is obsolete and needs to be rewritten */
            dp_ptr pcrumb = dp_mat_entry(crumb_mat, i, j);
            best = _get_max3_checked(cell->m+0, cell->m+1, cell->m+2, rank);
            if (_fmpz_vec_equal(best->vec, cell->m[0].vec, rank))
            {
                *pcrumb |= CRUMB_TOP;
            }
            if (_fmpz_vec_equal(best->vec, cell->m[1].vec, rank))
            {
                *pcrumb |= CRUMB_DIAG;
            }
            if (_fmpz_vec_equal(best->vec, cell->m[2].vec, rank))
            {
                *pcrumb |= CRUMB_LEFT;
            }

            /*
            flint_printf("%wd %wd %wd %wd ", i, j, k, l);
            flint_printf("%g %g %g\n",
                    exp(cell->m0),
                    exp(cell->m1),
                    exp(cell->m2));
            */

            l += 2;
            i--;
            j++;
        }
    }

    /* report the probability matrices */
    /*
    slong c;
    arf_ptr mid;
    for (c = 0; c < 3; c++)
    {
        flint_printf("m%wd:\n", c);
        for (i = 0; i < nrows; i++)
        {
            for (j = 0; j < ncols; j++)
            {
                k = i + j;
                l = nrows - 1 + j - i;
                cell = hwave_mat_entry(wave, k, l);
                mid = arb_midref(cell->m[c].value);
                flint_printf("%.17lf ", exp(arf_get_d(mid, ARF_RND_NEAR)));
            }
            flint_printf("\n");
        }
        flint_printf("\n");
    }
    flint_printf("\n");
    */

    /* report the score */
    i = nrows - 1;
    j = ncols - 1;
    k = i + j;
    l = nrows - 1 + j - i;
    cell = hwave_mat_entry(wave, k, l);
    best = _get_max3_checked(cell->m+0, cell->m+1, cell->m+2, rank);
    _arb_vec_dot_fmpz_vec(sol->log_probability, v, best->vec, rank, prec);

    /* do the traceback */
    dp_mat_get_alignment(
            sol->A, sol->B, &(sol->len),
            crumb_mat, A, B);

    /* clear the tables */
    hwave_mat_clear(wave);
    dp_mat_clear(crumb_mat);
}


void
tkf91_dp_r(
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    /* Inputs:
     *   mat : the generator matrix -- mat_ij where i is a generator index
     *         and j is an expression index.
     *   expressions_table : map from expression index to expression object
     *   g : a struct with tkf91 generator indices
     */

    fmpz_mat_t H, V;
    slong level, prec, rank;
    tkf91_generator_vecs_t h;

    /*
     * For now, use an arbitrary precision.
     * In later versions of this program,
     * we could decide to adjust it dynamically if it is detected
     * to be insufficient.
     */
    level = 8;
    prec = 1 << level;

    /* Compute a Hermite decomposition of the generator matrix. */
    /* U*mat = H ; U^-1 = V ; rank = rank(H) */
    fmpz_mat_init(H, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
    fmpz_mat_init(V, fmpz_mat_nrows(mat), fmpz_mat_nrows(mat));
    _fmpz_mat_hnf_inverse_transform(H, V, &rank, mat);

    /*
    flint_fprintf(stderr, "matrix rank revealed by the ");
    flint_fprintf(stderr, "Hermite normal form: %wd\n", rank);
    */

    /*
     * 'vecify' a structure with tkf91 generators,
     * by pointing the member variables to copies of rows of V.
     */
    tkf91_generator_vecs_init(h, g, V, rank);

    arb_ptr v;
    v = _arb_vec_init(rank);
    compute_hlogy(v, H, expressions_table, rank, level);


    /*
     * Begin checking some invariants.
     * Require that the result of the calculation does not depend on the basis.
     *
     * Initialize an arbitrary precision matrix.
     * This will be part of the calculation G*log(y)
     * which gives the score of each generator.
     */
    /*
    arb_mat_t arbG;
    arb_mat_init(arbG, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
    arb_mat_set_fmpz_mat(arbG, mat);
    */

    /*
     * Compute the generator logs G*log(y).
     */
    /*
    arb_mat_t generator_logs;
    arb_mat_init(generator_logs, fmpz_mat_nrows(mat), 1);
    arb_mat_mul(generator_logs, arbG, expression_logs, prec);
    */

    /*
    flint_printf("log probability for each generator:\n");
    arb_mat_printd(generator_logs, 15);
    flint_printf("\n");
    */

    /*
     * Compute dot products which should be equivalent.
     */
    /*
    flint_printf("a presumably equivalent calculation in a different basis\n");
    for (i = 0; i < fmpz_mat_nrows(mat); i++)
    {
        _arb_vec_dot_fmpz_vec(x, v, V->rows[i], rank, prec);

        flint_printf("first %wd columns of row %wd of V : ", rank, i);
        flint_printf("\n");
        _fmpz_vec_print(V->rows[i], rank);
        flint_printf("dot product : ");

        arb_printd(x, 15);
        flint_printf("\n");
    }
    flint_printf("\n");
    */

    /*
    flint_printf("V:\n"); fmpz_mat_print_pretty(V); flint_printf("\n");
    flint_printf("H:\n"); fmpz_mat_print_pretty(H); flint_printf("\n");
    */

    /* clear temporary variables */
    fmpz_mat_clear(H);
    fmpz_mat_clear(V);
    /*
    arb_mat_clear(arbG);
    arb_mat_clear(generator_logs);
    */

    /* do the thing */
    tkf91_dynamic_programming_hermite(
            sol, req, h, v, rank, prec, A, szA, B, szB);

    /* this will have aborted if not optimal so set optimality to true */
    sol->optimality_flag = 1;

    /* clear the remaining variables */
    _arb_vec_clear(v, rank);
    tkf91_generator_vecs_clear(h);
}
