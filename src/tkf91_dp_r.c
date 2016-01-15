/*
 * tkf91 dynamic programming cell bounds
 * using a dense tableau with arb_t real balls.
 */

#include <time.h>

#include "mag.h"
#include "arb.h"
#include "arf.h"
#include "arb_mat.h"

#include "tkf91_dp.h"
#include "tkf91_dp_r.h"
#include "tkf91_dp_bound.h"
#include "dp.h"
#include "forward.h"
#include "printutil.h"
#include "unused.h"
#include "bound_mat.h"


typedef struct
{
    arb_t m1_00;
    arb_t m0_10;
    arb_struct m0_i0_incr[4];
    arb_t m2_01;
    arb_struct m2_0j_incr[4];
    arb_struct c0_incr[4];
    arb_struct c1_incr[16];
    arb_struct c2_incr[4];
} tkf91_values_struct;
typedef tkf91_values_struct tkf91_values_t[1];
typedef tkf91_values_struct * tkf91_values_ptr;


static void tkf91_values_init(tkf91_values_t h,
        const tkf91_generator_indices_t g,
        arb_ptr m);

static void tkf91_values_clear(tkf91_values_t h);

static void _bounds_init(tkf91_values_t v, slong level,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g);

static void _arb_mat_get_col(arb_ptr v, const arb_mat_t mat, slong j);
static void _arb_max(arb_t z, const arb_t x, const arb_t y);
static void _arb_init_set(arb_t z, const arb_t x);

void
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
_arb_init_set(arb_t z, const arb_t x)
{
    arb_init(z);
    arb_set(z, x);
}

void
_arb_max(arb_t z, const arb_t x, const arb_t y)
{
    if (arb_lt(x, y))
    {
        arb_set(z, y);
    }
    else if (arb_lt(y, x))
    {
        arb_set(z, x);
    }
    else
    {
        arf_max(arb_midref(z), arb_midref(x), arb_midref(y));
        mag_max(arb_radref(z), arb_radref(x), arb_radref(y));
    }
}


void
tkf91_values_init(
        tkf91_values_t h,
        const tkf91_generator_indices_t g,
        arb_ptr m)
{
    slong i, j;
    _arb_init_set(h->m1_00, m+g->m1_00);
    _arb_init_set(h->m0_10, m+g->m0_10);
    _arb_init_set(h->m2_01, m+g->m2_01);
    for (i = 0; i < 4; i++)
    {
        _arb_init_set(h->m0_i0_incr+i, m+g->m0_i0_incr[i]);
        _arb_init_set(h->m2_0j_incr+i, m+g->m2_0j_incr[i]);
        _arb_init_set(h->c0_incr+i, m+g->c0_incr[i]);
        for (j = 0; j < 4; j++)
        {
            _arb_init_set(h->c1_incr+i*4+j, m+g->c1_incr[i*4+j]);
        }
        _arb_init_set(h->c2_incr+i, m+g->c2_incr[i]);
    }
}

void
tkf91_values_clear(tkf91_values_t h)
{
    slong i, j;
    arb_clear(h->m1_00);
    arb_clear(h->m0_10);
    arb_clear(h->m2_01);
    for (i = 0; i < 4; i++)
    {
        arb_clear(h->m0_i0_incr+i);
        arb_clear(h->m2_0j_incr+i);
        arb_clear(h->c0_incr+i);
        for (j = 0; j < 4; j++)
        {
            arb_clear(h->c1_incr+i*4+j);
        }
        arb_clear(h->c2_incr+i);
    }
}



void
_bounds_init(tkf91_values_t h, slong level,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g)
{
    slong nr, nc, prec;
    arb_mat_t G, U, V;

    prec = 1 << level;

    /* count the generators and expressions respectively */
    nr = fmpz_mat_nrows(mat);
    nc = fmpz_mat_ncols(mat);

    /* initialize the arbitrary precision integer exponent matrix */
    arb_mat_init(G, nr, nc);
    arb_mat_set_fmpz_mat(G, mat);

    /* compute logs of expressions to the specified precision */
    {
        slong i;
        arb_t x;
        arb_init(x);
        arb_mat_init(U, nc, 1);
        for (i = 0; i < nc; i++)
        {
            expr_eval(x, expressions_table[i], level);
            arb_log(arb_mat_entry(U, i, 0), x, prec);
        }
        arb_clear(x);
    }

    /* compute logs of generators */
    arb_mat_init(V, nr, 1);
    arb_mat_mul(V, G, U, prec);

    /* 
     * Copy the column vector to a new vector
     * and use it to initialize the structure.
     */
    {
        arb_ptr v;
        v = _arb_vec_init(nr);
        _arb_mat_get_col(v, V, 0);
        tkf91_values_init(h, g, v);
        _arb_vec_clear(v, nr);
    }

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
    arb_struct max2;
    arb_struct max3;
} cell_struct;
typedef cell_struct cell_t[1];
typedef cell_struct * cell_ptr;

static void cell_init(cell_t x);
static void cell_clear(cell_t x);

void
cell_init(cell_t x)
{
    arb_init(&(x->max2));
    arb_init(&(x->max3));
}

void
cell_clear(cell_t x)
{
    arb_clear(&(x->max2));
    arb_clear(&(x->max3));
}




/* the tableau cell visitor sees this data */
typedef struct
{
    slong level;
    arb_t m0;
    arb_t m1;
    arb_t m2;
    tkf91_values_t h;
    const slong *A;
    const slong *B;
} utility_struct;
typedef utility_struct utility_t[1];
typedef utility_struct * utility_ptr;

static void utility_clear(utility_t p);
static void utility_init(utility_t p,
        slong level,
        fmpz_mat_t mat,
        expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, const slong *B);

void
utility_init(utility_t p,
        slong level,
        fmpz_mat_t mat,
        expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, const slong *B)
{
    _bounds_init(p->h, level, mat, expressions_table, g);
    arb_init(p->m0);
    arb_init(p->m1);
    arb_init(p->m2);
    p->level = level;
    p->A = A;
    p->B = B;
}

void
utility_clear(utility_t p)
{
    tkf91_values_clear(p->h);
    arb_clear(p->m0);
    arb_clear(p->m1);
    arb_clear(p->m2);
}



static void *_init(void *userdata, size_t num);
static void _clear(void *userdata, void *celldata, size_t num);
static int _visit(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left);
static int _visit_boundary(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left);
static int _visit_center(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left);


void *
_init(void *userdata, size_t num)
{
    UNUSED(userdata);
    cell_ptr p = malloc(num * sizeof(cell_struct));
    size_t i;
    for (i = 0; i < num; i++)
    {
        cell_init(p + i);
    }
    return p;
}


void
_clear(void *userdata, void *celldata, size_t num)
{
    UNUSED(userdata);
    cell_ptr p = celldata;
    size_t i;
    for (i = 0; i < num; i++)
    {
        cell_clear(p + i);
    }
    free(p);
}


int
_visit_boundary(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    utility_ptr p = userdata;
    tkf91_values_ptr h = p->h;
    slong prec = 1 << p->level;
    UNUSED(curr);
    UNUSED(diag);
    UNUSED(mat);

    arb_neg_inf(p->m0);
    arb_neg_inf(p->m1);
    arb_neg_inf(p->m2);

    if (i == 0 && j == 0)
    {
        arb_set(p->m1, h->m1_00);
    }
    else if (i == 1 && j == 0)
    {
        arb_set(p->m0, h->m0_10);
    }
    else if (i == 0 && j == 1)
    {
        arb_set(p->m2, h->m2_01);
    }
    else
    {
        if (i == 0)
        {
            slong ntb = p->B[j - 1];
            cell_ptr p2 = left;
            arb_add(p->m2, &(p2->max2), h->m2_0j_incr+ntb, prec);
        }
        else if (j == 0)
        {
            slong nta = p->A[i - 1];
            cell_ptr p0 = top;
            arb_add(p->m0, &(p0->max3), h->m0_i0_incr+nta, prec);
        }
    }
    return 0;
}


int
_visit_center(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    utility_ptr p = userdata;
    tkf91_values_ptr h = p->h;
    dp_t x = *dp_mat_entry(mat, i, j);
    slong prec = 1 << p->level;
    UNUSED(curr);

    slong nta = p->A[i - 1];
    slong ntb = p->B[j - 1];

    if (dp_m0_is_interesting(x))
    {
        cell_ptr p0 = top;
        arb_add(p->m0, &(p0->max3), h->c0_incr+nta, prec);
    }

    if (dp_m1_is_interesting(x))
    {
        cell_ptr p1 = diag;
        arb_add(p->m1, &(p1->max3), h->c1_incr+nta*4 + ntb, prec);
    }

    if (dp_m2_is_interesting(x))
    {
        cell_ptr p2 = left;
        arb_add(p->m2, &(p2->max2), h->c2_incr+ntb, prec);
    }

    return 0;
}


int
_visit(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    utility_ptr p = userdata;
    cell_ptr c = curr;
    dp_t *px = dp_mat_entry(mat, i, j);
    dp_t x = *px;

    /*
     * Update the upper and lower bounds of the subset of m0, m1, m2
     * that is interesting for this cell according to the tableau flags.
     */
    if (i < 1 || j < 1)
    {
        _visit_boundary(userdata, mat, i, j, curr, top, diag, left);
    }
    else
    {
        _visit_center(userdata, mat, i, j, curr, top, diag, left);
    }

    /* If max2 is interesting for this cell then update its bounds. */
    if (x & DP_MAX2)
    {
        arb_neg_inf(&(c->max2));
        if (x & DP_MAX2_M1)
        {
            _arb_max(&(c->max2), &(c->max2), p->m1);
        }
        if (x & DP_MAX2_M2)
        {
            _arb_max(&(c->max2), &(c->max2), p->m2);
        }
    }

    /* If max3 is interesting for this cell then update its bounds. */
    if (x & DP_MAX3)
    {
        arb_neg_inf(&(c->max3));
        if (x & DP_MAX3_M0)
        {
            _arb_max(&(c->max3), &(c->max3), p->m0);
        }
        if (x & DP_MAX3_M1)
        {
            _arb_max(&(c->max3), &(c->max3), p->m1);
        }
        if (x & DP_MAX3_M2)
        {
            _arb_max(&(c->max3), &(c->max3), p->m2);
        }
    }

    /* If max2 is interesting for this cell then update its candidate flags */
    if (x & DP_MAX2)
    {
        if ((x & DP_MAX2_M1) && arb_lt(p->m1, &(c->max2)))
        {
            *px &= ~DP_MAX2_M1;
        }
        if ((x & DP_MAX2_M2) && arb_lt(p->m2, &(c->max2)))
        {
            *px &= ~DP_MAX2_M2;
        }
    }

    /* If max3 is interesting for this cell then update its candidate flags */
    if (x & DP_MAX3)
    {
        if ((x & DP_MAX3_M0) && arb_lt(p->m0, &(c->max3)))
        {
            *px &= ~DP_MAX3_M0;
        }
        if ((x & DP_MAX3_M1) && arb_lt(p->m1, &(c->max3)))
        {
            *px &= ~DP_MAX3_M1;
        }
        if ((x & DP_MAX3_M2) && arb_lt(p->m2, &(c->max3)))
        {
            *px &= ~DP_MAX3_M2;
        }
    }

    return 0;
}


void
tkf91_dp_r(
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    slong level = 8;
    tkf91_dp_r_level(level,
            sol, req, mat, expressions_table, g, A, szA, B, szB);
}


void
tkf91_dp_r_level(slong level,
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    utility_t util;
    forward_strategy_t s;
    slong nrows, ncols;
    clock_t start;
    int verbose = 0;
    FILE *file = NULL;
    if (verbose)
    {
        file = stderr;
    }

    if (!req->trace)
    {
        flint_fprintf(stderr, "tkf91_dp_r_(level %wd): ", level);
        flint_fprintf(stderr, "req->trace is required\n");
        abort();
    }

    if (!sol->mat)
    {
        flint_fprintf(stderr, "tkf91_dp_r_(level %wd): ", level);
        flint_fprintf(stderr, "sol->mat is required\n");
        abort();
    }

    nrows = dp_mat_nrows(sol->mat);
    ncols = dp_mat_ncols(sol->mat);
    if (nrows != (slong) szA + 1 ||
        ncols != (slong) szB + 1)
    {
        flint_fprintf(stderr, "tkf91_dp_r_(level %wd): ", level);
        flint_fprintf(stderr, "the sequence lengths are ");
        flint_fprintf(stderr, "incompatible with the tableau dimensions\n");
        abort();
    }

    start = clock();
    utility_init(util, level, mat, expressions_table, g, A, B);
    s->init = _init;
    s->clear = _clear;
    s->visit = _visit;
    s->sz_celldata = sizeof(cell_struct);
    s->userdata = util;
    dp_forward(sol->mat, s);
    utility_clear(util);
    _fprint_elapsed(file, "dynamic programming", clock() - start);

    /* update flags using a backward pass through the tableau */
    start = clock();
    dp_mat_backward(sol->mat);
    _fprint_elapsed(file, "backward algorithm pass", clock() - start);

    /* extract the alignment */
    start = clock();
    dp_mat_get_alignment(
            sol->A, sol->B, &(sol->len),
            sol->mat, A, B);
    _fprint_elapsed(file, "alignment traceback", clock() - start);
}


void
tkf91_dp_high(
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    slong level = -1;
    sol->optimality_flag = 0;
    while (!sol->optimality_flag)
    {
        if (level < 0)
        {
            tkf91_dp_mag(
                    sol, req, mat, expressions_table, g,
                    A, szA, B, szB);
            level = 6;
        }
        else
        {
            tkf91_dp_r_level(level,
                    sol, req, mat, expressions_table, g,
                    A, szA, B, szB);
            level++;
        }
        tkf91_dp_verify_symbolically(
                &sol->optimality_flag, 
                mat, g, sol->mat,
                expressions_table,
                A, B);
    }
}
