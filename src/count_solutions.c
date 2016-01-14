#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"

#include "count_solutions.h"
#include "forward.h"
#include "dp.h"
#include "unused.h"


static void *_init(void *userdata, size_t num);
static void _clear(void *userdata, void *celldata, size_t num);
static int _visit(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left);


void *_init(void *userdata, size_t num)
{
    UNUSED(userdata);
    return _fmpz_vec_init(num);
}


void _clear(void *userdata, void *celldata, size_t num)
{
    UNUSED(userdata);
    _fmpz_vec_clear((fmpz *) celldata, num);
}


int _visit(void *userdata, dp_mat_t mat,
        slong i, slong j,
        void *curr, void *top, void *diag, void *left)
{
    fmpz *p = curr;
    dp_t x = *dp_mat_entry(mat, i, j);
    fmpz_zero(p);
    if (x & DP_TRACE)
    {
        if (i == 0 && j == 0)
        {
            fmpz_one(p);
        }
        if (top && (x & DP_MAX3_M0))
        {
            fmpz_add(p, p, (fmpz *) top);
        }
        if (diag && (x & DP_MAX3_M1))
        {
            fmpz_add(p, p, (fmpz *) diag);
        }
        if (left && (x & DP_MAX3_M2))
        {
            fmpz_add(p, p, (fmpz *) left);
        }
    }

    /* special-case the visit to the bottom right corner cell */
    {
        slong nrows = dp_mat_nrows(mat);
        slong ncols = dp_mat_ncols(mat);
        if (i == nrows-1 && j == ncols-1)
        {
            fmpz_set((fmpz *) userdata, p);
        }
    }

    return 0;
}


void
count_solutions(fmpz_t res, dp_mat_t mat)
{
    forward_strategy_t s;

    s->init = _init;
    s->clear = _clear;
    s->visit = _visit;
    s->sz_celldata = sizeof(fmpz);
    s->userdata = res;

    dp_forward(mat, s);
}
