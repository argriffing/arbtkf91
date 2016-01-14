#include "flint/flint.h"

#include "tkf91_generator_vecs.h"
#include "forward.h"
#include "dp.h"
#include "bound_mat.h"


static void *_init(void *userdata, size_t num);
static void _clear(void *userdata, void *celldata, size_t num);
static int _visit(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left);

static fmpz * _pmax2(fmpz * celldata, slong rank);
static fmpz * _pmax3(fmpz * celldata, slong rank);
static void _cell_set_max2(fmpz *c, const fmpz *v, slong rank);
static void _cell_set_max3(fmpz *c, const fmpz *v, slong rank);
static int _visit_boundary(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
static int _visit_center(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)




/* the tableau cell visitor sees this data */
typedef struct
{
    tkf91_generator_vecs_ptr h;
    fmpz *m0;
    fmpz *m1;
    fmpz *m2;
    slong *A;
    slong *B;
} utility_struct;
typedef utility_struct utility_t[1];
typedef utility_struct * utility_ptr;

static void utility_clear(utility_t p);
static void utility_init(utility_t p,
        tkf91_generator_vecs_t h, slong *A, slong *B);

void
utility_init(utility_t p, tkf91_generator_vecs_t h, slong *A, slong *B)
{
    p->h = h;
    slong rank = tkf91_generator_vecs_rank(p->h);
    p->m0 = _fmpz_vec_init(rank);
    p->m1 = _fmpz_vec_init(rank);
    p->m2 = _fmpz_vec_init(rank);
    p->A = A;
    p->B = B;
}

void
utility_clear(utility_t p)
{
    slong rank = tkf91_generator_vecs_rank(p->h);
    _fmpz_vec_clear(p->m0, rank);
    _fmpz_vec_clear(p->m1, rank);
    _fmpz_vec_clear(p->m2, rank);
}









fmpz *
_pmax2(fmpz * celldata, slong rank)
{
    return celldata + 0 * rank;
}


fmpz *
_pmax3(fmpz * celldata, slong rank)
{
    return celldata + 1 * rank;
}

void
_cell_set_max2(fmpz *c, const fmpz *v, slong rank)
{
    _fmpz_vec_set(_pmax2(c, rank), v, rank);
}

void
_cell_set_max3(fmpz *c, const fmpz *v, slong rank)
{
    _fmpz_vec_set(_pmax3(c, rank), v, rank);
}



void *
_init(void *userdata, size_t num)
{
    /*
     * The userdata provides the rank of the generator matrix.
     * The 'num' argument is the number of buffer cells requested
     * by the general forward algorithm framework.
     */
    utility_ptr p = userdata;
    slong rank = tkf91_generator_vecs_rank(p->h);
    slong len = (slong) (2 * rank * num);
    return _fmpz_vec_init(len);
}

void
_clear(void *userdata, void *celldata, size_t num)
{
    utility_ptr p = userdata;
    slong rank = tkf91_generator_vecs_rank(p->h);
    slong len = (slong) (2 * rank * num);
    _fmpz_vec_clear(celldata, len);
}

int
_visit_boundary(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    utility_ptr p = userdata;
    tkf91_generator_vecs_ptr h = p->h;
    slong rank = tkf91_generator_vecs_rank(h);

    _fmpz_vec_zero(p->m0, rank);
    _fmpz_vec_zero(p->m1, rank);
    _fmpz_vec_zero(p->m2, rank);
    if (i == 0 && j == 0)
    {
        _fmpz_vec_set(p->m1, h->m1_00, rank);
    }
    else if (i == 1 && j == 0)
    {
        _fmpz_vec_set(p->m0, h->m0_10, rank);
    }
    else if (i == 0 && j == 1)
    {
        _fmpz_vec_set(p->m2, h->m2_01, rank);
    }
    else
    {
        if (i == 0)
        {
            slong ntb = p->B[j - 1];
            cell_ptr p2 = node->left;
            fmpz *p2_max2 = _pmax2(p2, rank);
            _fmpz_vec_add(p->m2, p2_max2, h->m2_0j_incr[ntb], rank);
        }
        else if (j == 0)
        {
            slong nta = p->A[i - 1];
            cell_ptr p0 = node->top;
            fmpz *p0_max3 = _pmax3(p0, rank);
            _fmpz_vec_add(p->m0, p0_max3, h->m0_i0_incr[nta], rank);
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
    tkf91_generator_vecs_ptr h = p->h;
    slong rank = tkf91_generator_vecs_rank(h);
    dp_t x = *dp_mat_entry(mat, i, j);

    slong nta = p->A[i - 1];
    slong ntb = p->B[j - 1];

    if (dp_m0_is_interesting(x))
    {
        fmpz *p0_max3 = _pmax3(node->top, rank);
        _fmpz_vec_add(p->m0, p0_max3, h->c0_incr[nta], rank);
    }

    if (dp_m1_is_interesting(x))
    {
        fmpz *p1_max3 = _pmax3(node->diag, rank);
        _fmpz_vec_add(p->m1, p1_max3, h->c1_incr[nta*4+ntb], rank);
    }

    if (dp_m2_is_interesting(x))
    {
        fmpz *p2_max2 = _pmax2(node->left, rank);
        _fmpz_vec_add(p->m2, p2_max2, h->c2_incr[ntb], rank);
    }
    return 0;
}


int
_visit_check_consensus(
        void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    utility_ptr p = userdata;
    tkf91_generator_vecs_ptr h = p->h;
    slong rank = tkf91_generator_vecs_rank(h);
    dp_t x = *dp_mat_entry(mat, i, j);

    if (x & DP_MAX2)
    {
        if ((x & DP_MAX2_M1) && (x & DP_MAX2_M2))
        {
            if (!_fmpz_vec_equal(p->m1, p->m2, rank)) return -1;
        }
    }
    if (x & DP_MAX3)
    {
        if ((x & DP_MAX3_M0) && (x & DP_MAX3_M1))
        {
            if (!_fmpz_vec_equal(p->m0, p->m1, rank)) return -1;
        }
        if ((x & DP_MAX3_M1) && (x & DP_MAX3_M2))
        {
            if (!_fmpz_vec_equal(p->m1, p->m2, rank)) return -1;
        }
        if ((x & DP_MAX3_M2) && (x & DP_MAX3_M0))
        {
            if (!_fmpz_vec_equal(p->m2, p->m0, rank)) return -1;
        }
    }
    return 0;
}



int
_visit_update_celldata(
        void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    utility_ptr p = userdata;
    tkf91_generator_vecs_ptr h = p->h;
    slong rank = tkf91_generator_vecs_rank(h);
    dp_t x = *dp_mat_entry(mat, i, j);

    if (x & DP_MAX2)
    {
        if (x & DP_MAX2_M1)
        {
            _cell_set_max2(curr, p->m1, rank);
        }
        else if (x & DP_MAX2_M2)
        {
            _cell_set_max2(curr, p->m2, rank);
        }
    }
    if (x & DP_MAX3)
    {
        if (x & DP_MAX3_M0)
        {
            _cell_set_max3(curr, p->m0, rank);
        }
        else if (x & DP_MAX3_M1)
        {
            _cell_set_max3(curr, p->m1, rank);
        }
        else if (x & DP_MAX3_M2)
        {
            _cell_set_max3(curr, p->m2, rank);
        }
    }
    return 0;
}


int
_visit(
        void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    int result;

    /*
     * Compute vectors associated with the subset of {m0, m1, m2}
     * that is considered to be interesting for this cell according
     * to the tableau flags.
     * These vector values are stored in the object referenced
     * by the utility pointer.
     */
    if (i < 1 || j < 1)
    {
        _visit_boundary(userdata, mat, i, j, curr, top, diag, left);
    }
    else
    {
        _visit_center(userdata, mat, i, j, curr, top, diag, left);
    }

    /*
     * For max2 and max3, if the max is interesting then check
     * for consensus among candidates.
     * If there is no consensus then return a nonzero integer.
     */
    result = _visit_check_consensus(userdata, mat, i, j, curr, top, diag, left);
    if (result)
    {
        return result;
    }

    /* Update the cell data. */
    _visit_update_celldata(userdata, mat, i, j, curr, top, diag, left);

    return 0;
}







int
tkf91_dp_verify_symbolically(
        int *verified,
        fmpz_mat_t mat,
        const tkf91_generator_indices_t g,
        dp_mat_t tableau,
        expr_ptr * expressions_table,
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
    arb_ptr v;
    slong level;

    level = 8;

    /* Compute a Hermite decomposition of the generator matrix. */
    /* U*mat = H ; U^-1 = V ; rank = rank(H) */
    fmpz_mat_init(H, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
    fmpz_mat_init(V, fmpz_mat_nrows(mat), fmpz_mat_nrows(mat));
    _fmpz_mat_hnf_inverse_transform(H, V, &rank, mat);

    tkf91_generator_vecs_init(h, g, V, rank);
    bound_mat_init(b, mask, rank);

    v = _arb_vec_init(rank);
    compute_hlogy(v, H, expressions_table, rank, level);

    /* this block replaces the old call to the symbolic verification */
    int result;
    {
        forward_strategy_t s;
        utility_t util;

        utility_init(util, h, A, B);
        s->init = _init;
        s->clear = _clear;
        s->visit = _visit;
        s->sz_celldata = sizeof(cell_struct);
        s->userdata = util;

        result = dp_forward(tableau, s);
        utility_clear(util);
    }

    _arb_vec_clear(v, rank);
    fmpz_mat_clear(H);
    fmpz_mat_clear(V);
    tkf91_generator_vecs_clear(h);

    return result;
}




/* the following functions are (were?) for debugging */

int _check_equal(fmpz * a, fmpz * b, slong r, arb_ptr v);

int
_check_equal(fmpz * a, fmpz * b, slong r, arb_ptr v)
{
    if (!_fmpz_vec_equal(a, b, r))
    {
        flint_printf("expected two vectors to be equal, but they aren't\n");

        arb_t x, y;
        arb_init(x);
        arb_init(y);
        _arb_vec_dot_fmpz_vec(x, v, a, r, 1 << 8);
        _arb_vec_dot_fmpz_vec(y, v, b, r, 1 << 8);

        _fmpz_vec_print(a, r);
        flint_printf(" : ");
        arb_print(x);
        flint_printf("\n");

        _fmpz_vec_print(b, r);
        flint_printf(" : ");
        arb_print(y);
        flint_printf("\n");

        flint_printf("overlap? %d\n", arb_overlaps(x, v));

        return 0;
    }
    return 1;
}

void _report(const char * name, fmpz * a, slong r);

void
_report(const char * name, fmpz * a, slong r)
{
    flint_printf(name);
    flint_printf(" integer vector : ");
    if (a)
    {
        _fmpz_vec_print(a, r);
    }
    else
    {
        flint_printf("unavailable");
    }
    flint_printf("\n");
}
