/*
 * tkf91 dynamic programming cell bounds
 * using a dense tableau with mag_t lower and upper bounds.
 */

#include <time.h>

#include "mag.h"
#include "arb_mat.h"

#include "tkf91_dp.h"
#include "tkf91_dp_bound.h"
#include "breadcrumbs.h"
#include "bound_mat.h"
#include "printutil.h"
#include "vis.h"
#include "count_solutions.h"


typedef struct
{
    mag_t m1_00;
    mag_t m0_10;
    mag_struct m0_i0_incr[4];
    mag_t m2_01;
    mag_struct m2_0j_incr[4];
    mag_struct c0_incr[4];
    mag_struct c1_incr[16];
    mag_struct c2_incr[4];
} tkf91_values_struct;
typedef tkf91_values_struct tkf91_values_t[1];


void tkf91_values_init(tkf91_values_t h,
        const tkf91_generator_indices_t g,
        mag_ptr m);

void tkf91_values_clear(tkf91_values_t h);

void _bounds_init(tkf91_values_t lb, tkf91_values_t ub,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g);


static __inline__ void
_arb_mat_get_col(arb_ptr v, const arb_mat_t mat, slong j)
{
    slong i, nr;
    nr = arb_mat_nrows(mat);
    for (i = 0; i < nr; i++)
    {
        arb_set(v+i, arb_mat_entry(mat, i, j));
    }
}



void
tkf91_values_init(
        tkf91_values_t h,
        const tkf91_generator_indices_t g,
        mag_ptr m)
{
    slong i, j;
    mag_init_set(h->m1_00, m+g->m1_00);
    mag_init_set(h->m0_10, m+g->m0_10);
    mag_init_set(h->m2_01, m+g->m2_01);
    for (i = 0; i < 4; i++)
    {
        mag_init_set(h->m0_i0_incr+i, m+g->m0_i0_incr[i]);
        mag_init_set(h->m2_0j_incr+i, m+g->m2_0j_incr[i]);
        mag_init_set(h->c0_incr+i, m+g->c0_incr[i]);
        for (j = 0; j < 4; j++)
        {
            mag_init_set(h->c1_incr+i*4+j, m+g->c1_incr[i*4+j]);
        }
        mag_init_set(h->c2_incr+i, m+g->c2_incr[i]);
    }
}

void
tkf91_values_clear(tkf91_values_t h)
{
    slong i, j;
    mag_clear(h->m1_00);
    mag_clear(h->m0_10);
    mag_clear(h->m2_01);
    for (i = 0; i < 4; i++)
    {
        mag_clear(h->m0_i0_incr+i);
        mag_clear(h->m2_0j_incr+i);
        mag_clear(h->c0_incr+i);
        for (j = 0; j < 4; j++)
        {
            mag_clear(h->c1_incr+i*4+j);
        }
        mag_clear(h->c2_incr+i);
    }
}



void
_bounds_init(tkf91_values_t lb, tkf91_values_t ub,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g)
{
    slong nr, nc, level, prec;
    arb_mat_t G, U, V;
    arb_t x;
    arb_ptr v;
    mag_ptr lb_arr, ub_arr;
    slong i;

    /* debug change this precision level */

    /* set the precision level used before converting to 30 bit magnitudes */
    level = 8;
    prec = 1 << level;

    /* count the generators and expressions respectively */
    nr = fmpz_mat_nrows(mat);
    nc = fmpz_mat_ncols(mat);

    /* initialize the arbitrary precision integer exponent matrix */
    arb_mat_init(G, nr, nc);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute logs of expressions to the specified precision */
    arb_init(x);
    arb_mat_init(U, nc, 1);
    for (i = 0; i < nc; i++)
    {
        expr_eval(x, expressions_table[i], level);
        arb_log(arb_mat_entry(U, i, 0), x, prec);
    }
    arb_clear(x);

    /* compute logs of generators */
    arb_mat_init(V, nr, 1);
    arb_mat_mul(V, G, U, prec);

    /* copy the column vector to a new vector and exponentiate its entries */
    v = _arb_vec_init(nr);
    _arb_mat_get_col(v, V, 0);
    for (i = 0; i < nr; i++)
    {
        arb_exp(v+i, v+i, prec);
    }

    /* create the lb and ub arrays */
    lb_arr = _mag_vec_init(nr);
    ub_arr = _mag_vec_init(nr);
    for (i = 0; i < nr; i++)
    {
        arb_get_mag_lower(lb_arr+i, v+i);
        arb_get_mag(ub_arr+i, v+i);
    }

    /* initialize the lb and ub structures */
    tkf91_values_init(lb, g, lb_arr);
    tkf91_values_init(ub, g, ub_arr);

    _mag_vec_clear(lb_arr, nr);
    _mag_vec_clear(ub_arr, nr);
    _arb_vec_clear(v, nr);
    arb_mat_clear(G);
    arb_mat_clear(U);
    arb_mat_clear(V);
}


/*
 * Each cell stores the lower bound and upper bound
 * for max(m1, m2) and max(m0, m1, m2).
 */

typedef struct
{
    mag_struct lb2;
    mag_struct ub2;
    mag_struct lb3;
    mag_struct ub3;
} cell_struct;
typedef cell_struct cell_t[1];
typedef cell_struct * cell_ptr;

void cell_init(cell_t x);
void cell_clear(cell_t x);
void cell_get_crumb(breadcrumb_t * pcrumb, const cell_t x,
        const mag_t ub_m0, const mag_t ub_m1, const mag_t ub_m2,
        slong i, slong j);

static __inline__ void
cell_fill(cell_t x,
        const mag_t lb_m0, const mag_t lb_m1, const mag_t lb_m2,
        const mag_t ub_m0, const mag_t ub_m1, const mag_t ub_m2)
{
    mag_max(&(x->ub2), ub_m1, ub_m2);
    mag_max(&(x->lb2), lb_m1, lb_m2);
    mag_max(&(x->ub3), &(x->ub2), ub_m0);
    mag_max(&(x->lb3), &(x->lb2), lb_m0);
}

void
cell_init(cell_t x)
{
    mag_init(&(x->lb2));
    mag_init(&(x->ub2));
    mag_init(&(x->lb3));
    mag_init(&(x->ub3));
}

void
cell_clear(cell_t x)
{
    mag_clear(&(x->lb2));
    mag_clear(&(x->ub2));
    mag_clear(&(x->lb3));
    mag_clear(&(x->ub3));
}

void
cell_get_crumb(breadcrumb_t * pcrumb, const cell_t x,
        const mag_t ub_m0, const mag_t ub_m1, const mag_t ub_m2,
        slong i, slong j)
{
    *pcrumb = 0;

    /* breadcrumb for alignment traceback */
    if (i || j)
    {
        if (mag_cmp(&(x->lb3), ub_m0) <= 0)
        {
            *pcrumb |= CRUMB_TOP;
        }
        if (mag_cmp(&(x->lb3), ub_m1) <= 0)
        {
            *pcrumb |= CRUMB_DIAG;
        }
        if (mag_cmp(&(x->lb3), ub_m2) <= 0)
        {
            *pcrumb |= CRUMB_LEFT;
        }
    }

    /* breadcrumb for info propagation for symbolic optimality verification */
    if (j)
    {
        if (mag_cmp(&(x->lb2), ub_m1) <= 0)
        {
            *pcrumb |= CRUMB_DIAG2;
        }
        if (mag_cmp(&(x->lb2), ub_m2) <= 0)
        {
            *pcrumb |= CRUMB_LEFT2;
        }
    }
}


typedef struct
{
    cell_ptr data;
    slong r;
    slong c;
} cellfront_struct;
typedef cellfront_struct cellfront_t[1];
typedef cellfront_struct * cellfront_ptr;

void cellfront_init(cellfront_t x, slong nrows, slong ncols);
void cellfront_clear(cellfront_t x);

static __inline__ cell_ptr
cellfront_entry(cellfront_t x, slong i, slong j)
{
    return x->data + (i % 2) * (x->c) + j;
}

static __inline__ cell_ptr
cellfront_entry_top(cellfront_t x, slong i, slong j)
{
    return cellfront_entry(x, i-1, j);
}

static __inline__ cell_ptr
cellfront_entry_diag(cellfront_t x, slong i, slong j)
{
    return cellfront_entry(x, i-1, j-1);
}

static __inline__ cell_ptr
cellfront_entry_left(cellfront_t x, slong i, slong j)
{
    return cellfront_entry(x, i, j-1);
}


void
cellfront_init(cellfront_t x, slong nrows, slong ncols)
{
    slong i;
    x->r = nrows;
    x->c = ncols;
    x->data = flint_malloc(2 * ncols * sizeof(cell_struct));
    for (i = 0; i < 2 * ncols; i++)
    {
        cell_init(x->data+i);
    }
}

void
cellfront_clear(cellfront_t x)
{
    slong i;
    for (i = 0; i < 2 * x->c; i++)
    {
        cell_clear(x->data+i);
    }
    flint_free(x->data);
}



void
tkf91_dp_bound(
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    clock_t start;
    int verbose = 0;
    FILE * file = NULL;
    if (verbose)
    {
        file = stderr;
    }

    /* Convert the first few args to generator lower and upper bounds. */
    tkf91_values_t lb;
    tkf91_values_t ub;
    _bounds_init(lb, ub, mat, expressions_table, g);

    /*
     * Define the matrix to be used for the traceback.
     * The number of rows is one greater than the length of the first sequence,
     * and the number of columns is one greater than the length of the second
     * sequence.
     */
    slong nrows = szA + 1;
    slong ncols = szB + 1;
    breadcrumb_mat_t crumb_mat;
    breadcrumb_ptr pcrumb;

    /*
     * TODO for efficiency do not duplicate
     * crumb_mat and sol->pmask.
     */

    if (req->trace)
    {
        breadcrumb_mat_init(crumb_mat, nrows, ncols);
    }

    /* define the cellfront matrix */
    cellfront_t cells;
    cellfront_init(cells, nrows, ncols);

    start = clock();

    /* iterate over rows of the dynamic programming matrix */
    cell_ptr cell, p0, p1, p2;
    slong i, j;
    slong nta, ntb;
    mag_t lb_m0, lb_m1, lb_m2;
    mag_t ub_m0, ub_m1, ub_m2;
    mag_init(lb_m0); mag_init(lb_m1); mag_init(lb_m2);
    mag_init(ub_m0); mag_init(ub_m1); mag_init(ub_m2);
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            cell = cellfront_entry(cells, i, j);
            if (i < 1 || j < 1)
            {
                mag_zero(lb_m0);
                mag_zero(lb_m1);
                mag_zero(lb_m2);
                mag_zero(ub_m0);
                mag_zero(ub_m1);
                mag_zero(ub_m2);
                if (i == 0 && j == 0)
                {
                    mag_set(lb_m1, lb->m1_00);
                    mag_set(ub_m1, ub->m1_00);
                }
                else if (i == 1 && j == 0)
                {
                    mag_set(lb_m0, lb->m0_10);
                    mag_set(ub_m0, ub->m0_10);
                }
                else if (i == 0 && j == 1)
                {
                    mag_set(lb_m2, lb->m2_01);
                    mag_set(ub_m2, ub->m2_01);
                }
                else
                {
                    if (i == 0)
                    {
                        ntb = B[j - 1];
                        p2 = cellfront_entry_left(cells, i, j);
                        mag_mul_lower(lb_m2, &(p2->lb2), lb->m2_0j_incr+ntb);
                        mag_mul(ub_m2, &(p2->ub2), ub->m2_0j_incr+ntb);
                    }
                    else if (j == 0)
                    {
                        nta = A[i - 1];
                        p0 = cellfront_entry_top(cells, i, j);
                        mag_mul_lower(lb_m0, &(p0->lb3), lb->m0_i0_incr+nta);
                        mag_mul(ub_m0, &(p0->ub3), ub->m0_i0_incr+nta);
                    }
                }
            }
            else
            {
                nta = A[i - 1];
                ntb = B[j - 1];
                p0 = cellfront_entry_top(cells, i, j);
                p1 = cellfront_entry_diag(cells, i, j);
                p2 = cellfront_entry_left(cells, i, j);

                mag_mul_lower(lb_m0, &(p0->lb3), lb->c0_incr+nta);
                mag_mul(ub_m0, &(p0->ub3), ub->c0_incr+nta);

                mag_mul_lower(lb_m1, &(p1->lb3), lb->c1_incr+nta*4 + ntb);
                mag_mul(ub_m1, &(p1->ub3), ub->c1_incr+nta*4 + ntb);

                mag_mul_lower(lb_m2, &(p2->lb2), lb->c2_incr+ntb);
                mag_mul(ub_m2, &(p2->ub2), ub->c2_incr+ntb);
            }

            /* fill the current cell using the current m* bounds */
            cell_fill(cell, lb_m0, lb_m1, lb_m2, ub_m0, ub_m1, ub_m2);

            /* optionally fill the table for traceback */
            if (req->trace)
            {
                pcrumb = breadcrumb_mat_entry(crumb_mat, i, j);

                /* pass i and j to detect boundaries only */
                cell_get_crumb(pcrumb, cell, ub_m0, ub_m1, ub_m2, i, j);
            }
        }
    }

    _fprint_elapsed(file, "dynamic programming", clock() - start);

    if (verbose)
    {
        if (req->trace)
        {
            flint_fprintf(stdout, "mask before tracing:\n");
            breadcrumb_mat_fprint(stdout, crumb_mat);
            flint_fprintf(stdout, "\n");
        }
        else
        {
            flint_fprintf(stdout, "not showing mask before tracing\n");
        }
    }

    /* report the score */

    cell = cellfront_entry(cells, nrows - 1, ncols - 1);

    if (verbose)
    {
        flint_printf("best alignment lower bound probability: ");
        mag_print(&(cell->lb3));
        flint_printf("\n");

        flint_printf("best alignment upper bound probability: ");
        mag_print(&(cell->ub3));
        flint_printf("\n");
    }

    if (req->trace)
    {
        /*
         * Get the not-ruled-out cells of the dynamic programming table,
         * and get the cells that are not ruled out for contributing
         * information towards the optimal alignment.
         */
        start = clock();
        breadcrumb_mat_get_mask(crumb_mat, crumb_mat);
        _fprint_elapsed(file, "mask extraction", clock() - start);

        /* create the tableau png image */
        if (req->png_filename)
        {
            start = clock();
            write_tableau_image(req->png_filename, crumb_mat, "tkf91 tableau");
            _fprint_elapsed(file, "create tableau png", clock() - start);
        }

        /* do the traceback */
        start = clock();

        breadcrumb_mat_get_alignment(
                sol->A, sol->B, &(sol->len),
                crumb_mat, A, B);

        _fprint_elapsed(file, "alignment traceback", clock() - start);

        start = clock();
        int verified = 0;
        tkf91_dp_verify_symbolically(&verified, mat, g, crumb_mat, A, B);
        _fprint_elapsed(file, "symbolic verification", clock() - start);
        sol->optimality_flag = verified;

        {
            fmpz_t solution_count;
            fmpz_init(solution_count);

            start = clock();
            count_solutions(solution_count, crumb_mat);
            _fprint_elapsed(file, "counting solutions", clock() - start);

            fmpz_set(sol->best_tie_count, solution_count);
            sol->has_best_tie_count = 1;

            fmpz_clear(solution_count);
        }
    }

    cellfront_clear(cells);

    if (req->trace && sol->pmask)
    {
        breadcrumb_mat_set(sol->pmask, crumb_mat);
    }

    if (req->trace)
    {
        breadcrumb_mat_clear(crumb_mat);
    }

    tkf91_values_clear(lb);
    tkf91_values_clear(ub);

    mag_clear(lb_m0); mag_clear(lb_m1); mag_clear(lb_m2);
    mag_clear(ub_m0); mag_clear(ub_m1); mag_clear(ub_m2);
}
