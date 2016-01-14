/*
 * tkf91 dynamic programming cell bounds
 * using a dense tableau with mag_t lower and upper bounds.
 */

#include <time.h>

#include "mag.h"
#include "arb_mat.h"

#include "tkf91_dp.h"
#include "tkf91_dp_bound.h"
#include "dp.h"
#include "forward.h"
#include "printutil.h"


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

static void cell_init(cell_t x);
static void cell_clear(cell_t x);

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




/* the tableau cell visitor sees this data */
typedef struct
{
    mag_t lb_m0;
    mag_t lb_m1;
    mag_t lb_m2;
    mag_t ub_m0;
    mag_t ub_m1;
    mag_t ub_m2;
    tkf91_values_t lb;
    tkf91_values_t ub;
    slong *A;
    slong *B;
} utility_struct;
typedef utility_struct utility_t[1];
typedef utility_struct * utility_ptr;

static void utility_clear(utility_t p);
static void utility_init(utility_t p,
        fmpz_mat_t mat,
        expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        slong *A, slong *B)

void
utility_init(utility_t p,
        fmpz_mat_t mat,
        expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        slong *A, slong *B)
{
    _bounds_init(p->lb, p->ub, mat, expressions_table, g);
    mag_init(p->lb_m0);
    mag_init(p->lb_m1);
    mag_init(p->lb_m2);
    mag_init(p->ub_m0);
    mag_init(p->ub_m1);
    mag_init(p->ub_m2);
    p->A = A;
    p->B = B;
}

void
utility_clear(utility_t p)
{
    tkf91_values_clear(p->lb);
    tkf91_values_clear(p->ub);
    mag_clear(p->lb_m0);
    mag_clear(p->lb_m1);
    mag_clear(p->lb_m2);
    mag_clear(p->ub_m0);
    mag_clear(p->ub_m1);
    mag_clear(p->ub_m2);
}



static void *_init(void *userdata, size_t num);
static void _clear(void *userdata, void *celldata, size_t num);
static int _visit(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left);
static int _visit_boundary(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
static int _visit_center(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)


void *
_init(void *userdata, size_t num)
{
    cell_ptr *p = malloc(num * sizeof(cell_struct));
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
    cell_ptr *p = celldata;
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

    mag_zero(p->lb_m0);
    mag_zero(p->lb_m1);
    mag_zero(p->lb_m2);
    mag_zero(p->ub_m0);
    mag_zero(p->ub_m1);
    mag_zero(p->ub_m2);

    if (i == 0 && j == 0)
    {
        mag_set(p->lb_m1, p->lb->m1_00);
        mag_set(p->ub_m1, p->ub->m1_00);
    }
    else if (i == 1 && j == 0)
    {
        mag_set(p->lb_m0, p->lb->m0_10);
        mag_set(p->ub_m0, p->ub->m0_10);
    }
    else if (i == 0 && j == 1)
    {
        mag_set(p->lb_m2, p->lb->m2_01);
        mag_set(p->ub_m2, p->ub->m2_01);
    }
    else
    {
        if (i == 0)
        {
            slong ntb = p->B[j - 1];
            cell_ptr p2 = left;
            mag_mul_lower(p->lb_m2, &(p2->lb2), p->lb->m2_0j_incr+ntb);
            mag_mul(p->ub_m2, &(p2->ub2), p->ub->m2_0j_incr+ntb);
        }
        else if (j == 0)
        {
            slong nta = p->A[i - 1];
            cell_ptr p0 = top;
            mag_mul_lower(p->lb_m0, &(p0->lb3), p->lb->m0_i0_incr+nta);
            mag_mul(p->ub_m0, &(p0->ub3), p->ub->m0_i0_incr+nta);
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
    dp_t x = *dp_mat_entry(mat, i, j);

    slong nta = p->A[i - 1];
    slong ntb = p->B[j - 1];

    if (dp_m0_is_interesting(x))
    {
        cell_ptr p0 = top;
        mag_mul_lower(p->lb_m0, &(p0->lb3), p->lb->c0_incr+nta);
        mag_mul(p->ub_m0, &(p0->ub3), p->ub->c0_incr+nta);
    }

    if (dp_m1_is_interesting(x))
    {
        cell_ptr p1 = diag;
        mag_mul_lower(p->lb_m1, &(p1->lb3), p->lb->c1_incr+nta*4 + ntb);
        mag_mul(p->ub_m1, &(p1->ub3), p->ub->c1_incr+nta*4 + ntb);
    }

    if (dp_m2_is_interesting(x))
    {
        cell_ptr p2 = left;
        mag_mul_lower(p->lb_m2, &(p2->lb2), p->lb->c2_incr+ntb);
        mag_mul(p->ub_m2, &(p2->ub2), p->ub->c2_incr+ntb);
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
        mag_zero(&(c->ub2));
        mag_zero(&(c->lb2));
        if (x & DP_MAX2_M1)
        {
            mag_max(&(c->ub2), &(c->ub2), p->ub_m1);
            mag_max(&(c->lb2), &(c->lb2), p->lb_m1);
        }
        if (x & DP_MAX2_M2)
        {
            mag_max(&(c->ub2), &(c->ub2), p->ub_m2);
            mag_max(&(c->lb2), &(c->lb2), p->lb_m2);
        }
    }

    /* If max3 is interesting for this cell then update its bounds. */
    if (x & DP_MAX3)
    {
        mag_zero(&(c->ub3));
        mag_zero(&(c->lb3));
        if (x & DP_MAX3_M0)
        {
            mag_max(&(c->ub3), &(c->ub3), p->ub_m0);
            mag_max(&(c->lb3), &(c->lb3), p->lb_m0);
        }
        if (x & DP_MAX3_M1)
        {
            mag_max(&(c->ub3), &(c->ub3), p->ub_m1);
            mag_max(&(c->lb3), &(c->lb3), p->lb_m1);
        }
        if (x & DP_MAX3_M2)
        {
            mag_max(&(c->ub3), &(c->ub3), p->ub_m2);
            mag_max(&(c->lb3), &(c->lb3), p->lb_m2);
        }
    }

    /* If max2 is interesting for this cell then update its candidate flags */
    if (x & DP_MAX2)
    {
        if ((x & DP_MAX2_M1) && (mag_cmp(p->ub_m1, &(c->lb2) < 0)))
        {
            *dp &= ~DP_MAX2_M1;
        }
        if ((x & DP_MAX2_M2) && (mag_cmp(p->ub_m2, &(c->lb2) < 0)))
        {
            *dp &= ~DP_MAX2_M2;
        }
    }

    /* If max3 is interesting for this cell then update its candidate flags */
    if (x & DP_MAX3)
    {
        if ((x & DP_MAX3_M0) && (mag_cmp(p->ub_m0, &(c->lb3) < 0)))
        {
            *dp &= ~DP_MAX3_M0;
        }
        if ((x & DP_MAX3_M1) && (mag_cmp(p->ub_m1, &(c->lb3) < 0)))
        {
            *dp &= ~DP_MAX3_M1;
        }
        if ((x & DP_MAX3_M2) && (mag_cmp(p->ub_m2, &(c->lb3) < 0)))
        {
            *dp &= ~DP_MAX3_M2;
        }
    }

    return 0;
}




void
tkf91_dp_bound(
        solution_t sol, const request_t req,
        fmpz_mat_t mat, expr_ptr * expressions_table,
        const tkf91_generator_indices_t g,
        const slong *A, size_t szA,
        const slong *B, size_t szB)
{
    utility_t util;
    forward_strategy_t s;
    clock_t start;
    int verbose = 0;
    FILE *file = NULL;
    if (verbose)
    {
        file = stderr;
    }

    if (!req->trace)
    {
        fprintf(stderr, "tkf91_dp_bound: req->trace is required\n");
        abort();
    }

    if (!sol->mat)
    {
        fprintf(stderr, "tkf91_dp_bound: sol->mat is required\n");
        abort();
    }

    start = clock();
    utility_init(util, mat, expressions_table, g, A, B);
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
