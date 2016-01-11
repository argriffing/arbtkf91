#include "flint/flint.h"
#include "flint/fmpq.h"
#include "arb.h"
#include "femtocas.h"


/* generic expression initialization */

void _expr_init(expr_ptr x);
void _default_clear(expr_data_ptr p);
void _default_print(expr_data_ptr p);
void _default_eval(arb_t res, expr_data_ptr p, slong level);


/* type-specific expression member functions */

void _expr_i_init(expr_ptr x, const fmpz_t a);
void _expr_i_clear(expr_data_ptr data);
void _expr_i_print(expr_data_ptr data);

void _expr_q_init(expr_ptr x, const fmpq_t a);
void _expr_q_clear(expr_data_ptr data);
void _expr_q_print(expr_data_ptr data);

void _expr_x_init(expr_ptr x, expr_ptr a);
void _expr_x_clear(expr_data_ptr data);
void _expr_x_print(expr_data_ptr data);

void _expr_xx_init(expr_ptr x, expr_ptr a, expr_ptr b);
void _expr_xx_clear(expr_data_ptr data);
void _expr_xx_print(expr_data_ptr data);


void
_default_clear(expr_data_ptr p)
{
    UNUSED(p);
    flint_printf("clear() not implemented\n");
    abort();
}

void
_default_print(expr_data_ptr p)
{
    UNUSED(p);
    flint_printf("print() not implemented\n");
    abort();
}

void
_default_eval(arb_t res, expr_data_ptr p, slong level)
{
    UNUSED(res);
    UNUSED(p);
    UNUSED(level);
    flint_printf("eval() not implemented\n");
    abort();
}

void
_expr_init(expr_ptr x)
{
    slong i;
    for (i = 0; i < EXPR_CACHE_CAP; i++)
    {
        arb_init(x->cache+i);
    }
    x->cachesize = 0;
    x->data = NULL;
    x->clear = &_default_clear;
    x->print = &_default_print;
    x->eval = &_default_eval;
}

void
expr_clear(expr_ptr x)
{
    slong i;
    for (i = 0; i < EXPR_CACHE_CAP; i++)
    {
        arb_clear(x->cache+i);
    }
    x->cachesize = 0;
    x->clear(x->data);
    x->data = NULL;
    x->clear = &_default_clear;
    x->print = &_default_print;
    x->eval = &_default_eval;
}

void
expr_print(expr_ptr x)
{
    x->print(x->data);
}

void
expr_eval(arb_t res, expr_ptr x, slong level)
{
    slong i;
    if (level < 0 || EXPR_CACHE_CAP <= level)
    {
        flint_printf("invalid log2 prec bits level %wd\n", level);
        abort();
    }
    if (level >= x->cachesize)
    {
        for (i = x->cachesize; i < level+1; i++)
        {
            x->eval(x->cache+i, x->data, i);
        }
        x->cachesize = level+1;
    }
    arb_set(res, x->cache+level);
}




/* implement signature-specific methods */

/*
 * typedef name codes:
 * i : fmpz
 * q : fmpq
 * x : expr_ptr
 *
 * Single letter codes like i, q, x indicate unary functions,
 * while double letter codes like xx indicate binary functions.
 */


/* i */

typedef struct
{
    fmpz a;
} _expr_i_struct;
typedef _expr_i_struct * _expr_i_ptr;

void
_expr_i_init(expr_ptr x, const fmpz_t a)
{
    _expr_i_ptr d = flint_malloc(sizeof(_expr_i_struct));
    fmpz_init_set(&(d->a), a);
    _expr_init(x);
    x->data = d;
    x->clear = &_expr_i_clear;
    x->print = &_expr_i_print;
}

void
_expr_i_clear(expr_data_ptr data)
{
    _expr_i_ptr d = (_expr_i_ptr) data;
    fmpz_clear(&(d->a));
    flint_free(d);
}

void
_expr_i_print(expr_data_ptr data)
{
    _expr_i_ptr d = (_expr_i_ptr) data;
    flint_printf("op("); fmpz_print(&(d->a)); flint_printf(")");
}





/* q */

typedef struct
{
    fmpq a;
} _expr_q_struct;
typedef _expr_q_struct * _expr_q_ptr;

void
_expr_q_init(expr_ptr x, const fmpq_t a)
{
    _expr_q_ptr d = flint_malloc(sizeof(_expr_q_struct));
    fmpq_init(&(d->a));
    fmpq_set(&(d->a), a);
    _expr_init(x);
    x->data = d;
    x->clear = &_expr_q_clear;
    x->print = &_expr_q_print;
}

void
_expr_q_clear(expr_data_ptr data)
{
    _expr_q_ptr d = (_expr_q_ptr) data;
    fmpq_clear(&(d->a));
    flint_free(d);
}

void
_expr_q_print(expr_data_ptr data)
{
    _expr_q_ptr d = (_expr_q_ptr) data;
    flint_printf("op("); fmpq_print(&(d->a)); flint_printf(")");
}





/* x */

typedef struct
{
    expr_ptr a;
} _expr_x_struct;
typedef _expr_x_struct * _expr_x_ptr;

void
_expr_x_init(expr_ptr x, expr_ptr a)
{
    _expr_x_ptr d = flint_malloc(sizeof(_expr_x_struct));
    d->a = a;
    _expr_init(x);
    x->data = d;
    x->clear = &_expr_x_clear;
    x->print = &_expr_x_print;
}

void
_expr_x_clear(expr_data_ptr data)
{
    _expr_x_ptr d = (_expr_x_ptr) data;
    d->a = NULL;
    flint_free(d);
}

void
_expr_x_print(expr_data_ptr data)
{
    _expr_x_ptr d = (_expr_x_ptr) data;
    flint_printf("op("); expr_print(d->a); flint_printf(")");
}





/* xx */

typedef struct
{
    expr_ptr a;
    expr_ptr b;
} _expr_xx_struct;
typedef _expr_xx_struct * _expr_xx_ptr;

void
_expr_xx_init(expr_ptr x, expr_ptr a, expr_ptr b)
{
    _expr_xx_ptr d = flint_malloc(sizeof(_expr_xx_struct));
    d->a = a;
    d->b = b;
    _expr_init(x);
    x->data = d;
    x->clear = &_expr_xx_clear;
    x->print = &_expr_xx_print;
}

void
_expr_xx_clear(expr_data_ptr data)
{
    _expr_xx_ptr d = (_expr_xx_ptr) data;
    d->a = NULL;
    d->b = NULL;
    flint_free(d);
}

void
_expr_xx_print(expr_data_ptr data)
{
    _expr_xx_ptr d = (_expr_xx_ptr) data;
    flint_printf("op("); expr_print(d->a);
    flint_printf(", "); expr_print(d->b); flint_printf(")");
}




/*
 * Constants and arithmetic are implemented here.
 * Each needs at least an evaluation function
 * and a constructor that points to that function.
 *
 * For most or all of these expressions,
 * the memory stuff is dealt with by the layer that is specific
 * to the type signature of the function.
 */



/* i */

void _expr_fmpz_eval(arb_t z, expr_data_ptr data, slong level);
void expr_fmpz(expr_ptr x, const fmpz_t a)
{
    _expr_i_init(x, a);
    x->eval = &_expr_fmpz_eval;
}
void _expr_fmpz_eval(arb_t z, expr_data_ptr data, slong level)
{
    UNUSED(level);
    _expr_i_ptr d = (_expr_i_ptr) data;
    arb_set_fmpz(z, &(d->a));
}

void _expr_log_fmpz_eval(arb_t z, expr_data_ptr data, slong level);
void expr_log_fmpz(expr_ptr x, const fmpz_t a)
{
    _expr_i_init(x, a);
    x->eval = &_expr_log_fmpz_eval;
}
void _expr_log_fmpz_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_i_ptr d = (_expr_i_ptr) data;
    arb_set_fmpz(z, &(d->a));
    arb_log(z, z, prec);
}


/* q */

void _expr_fmpq_eval(arb_t z, expr_data_ptr data, slong level);
void expr_fmpq(expr_ptr x, const fmpq_t a)
{
    _expr_q_init(x, a);
    x->eval = &_expr_fmpq_eval;
}
void _expr_fmpq_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_q_ptr d = (_expr_q_ptr) data;
    arb_set_fmpq(z, &(d->a), prec);
}

void _expr_log_fmpq_eval(arb_t z, expr_data_ptr data, slong level);
void expr_log_fmpq(expr_ptr x, const fmpq_t a)
{
    _expr_q_init(x, a);
    x->eval = &_expr_log_fmpq_eval;
}
void _expr_log_fmpq_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_q_ptr d = (_expr_q_ptr) data;
    arb_set_fmpq(z, &(d->a), prec);
    arb_log(z, z, prec);
}

void _expr_exp_fmpq_eval(arb_t z, expr_data_ptr data, slong level);
void expr_exp_fmpq(expr_ptr x, const fmpq_t a)
{
    _expr_q_init(x, a);
    x->eval = &_expr_exp_fmpq_eval;
}
void _expr_exp_fmpq_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_q_ptr d = (_expr_q_ptr) data;
    arb_set_fmpq(z, &(d->a), prec);
    arb_exp(z, z, prec);
}

void expr_complement_fmpq(expr_ptr x, const fmpq_t a)
{
    fmpq_t b;
    fmpq_init(b);
    fmpq_one(b);
    fmpq_sub(b, b, a);
    expr_fmpq(x, b);
    fmpq_clear(b);
}

void expr_neg_fmpq(expr_ptr x, const fmpq_t a)
{
    fmpq_t b;
    fmpq_init(b);
    fmpq_neg(b, a);
    expr_fmpq(x, b);
    fmpq_clear(b);
}


/* x */

void _expr_exp_eval(arb_t z, expr_data_ptr data, slong level);
void expr_exp(expr_ptr x, expr_ptr a)
{
    _expr_x_init(x, a);
    x->eval = &_expr_exp_eval;
}
void _expr_exp_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_x_ptr d = (_expr_x_ptr) data;
    expr_eval(z, d->a, level);
    arb_exp(z, z, prec);
}

void _expr_neg_eval(arb_t z, expr_data_ptr data, slong level);
void expr_neg(expr_ptr x, expr_ptr a)
{
    _expr_x_init(x, a);
    x->eval = &_expr_neg_eval;
}
void _expr_neg_eval(arb_t z, expr_data_ptr data, slong level)
{
    _expr_x_ptr d = (_expr_x_ptr) data;
    expr_eval(z, d->a, level);
    arb_neg(z, z);
}

void _expr_log_eval(arb_t z, expr_data_ptr data, slong level);
void expr_log(expr_ptr x, expr_ptr a)
{
    _expr_x_init(x, a);
    x->eval = &_expr_log_eval;
}
void _expr_log_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_x_ptr d = (_expr_x_ptr) data;
    expr_eval(z, d->a, level);
    arb_log(z, z, prec);
}

void _expr_log1p_eval(arb_t z, expr_data_ptr data, slong level);
void expr_log1p(expr_ptr x, expr_ptr a)
{
    _expr_x_init(x, a);
    x->eval = &_expr_log1p_eval;
}
void _expr_log1p_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_x_ptr d = (_expr_x_ptr) data;
    expr_eval(z, d->a, level);
    arb_log1p(z, z, prec);
}

void _expr_log1m_eval(arb_t z, expr_data_ptr data, slong level);
void expr_log1m(expr_ptr x, expr_ptr a)
{
    _expr_x_init(x, a);
    x->eval = &_expr_log1m_eval;
}
void _expr_log1m_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_x_ptr d = (_expr_x_ptr) data;
    expr_eval(z, d->a, level);
    arb_neg(z, z);
    arb_log1p(z, z, prec);
}

void _expr_complement_eval(arb_t z, expr_data_ptr data, slong level);
void expr_complement(expr_ptr x, expr_ptr a)
{
    _expr_x_init(x, a);
    x->eval = &_expr_complement_eval;
}
void _expr_complement_eval(arb_t z, expr_data_ptr data, slong level)
{
    slong prec = 1 << level;
    _expr_x_ptr d = (_expr_x_ptr) data;
    expr_eval(z, d->a, level);
    arb_t one;
    arb_init(one);
    arb_one(one);
    arb_sub(z, one, z, prec);
    arb_clear(one);
}


/* xx */


void _expr_add_eval(arb_t z, expr_data_ptr data, slong level);
void expr_add(expr_ptr x, expr_ptr a, expr_ptr b)
{
    _expr_xx_init(x, a, b);
    x->eval = &_expr_add_eval;
}
void _expr_add_eval(arb_t z, expr_data_ptr data, slong level)
{
    arb_t u, v;
    slong prec = 1 << level;
    _expr_xx_ptr d = (_expr_xx_ptr) data;
    arb_init(u);
    arb_init(v);
    expr_eval(u, d->a, level);
    expr_eval(v, d->b, level);
    arb_add(z, u, v, prec);
    arb_clear(u);
    arb_clear(v);
}

void _expr_mul_eval(arb_t z, expr_data_ptr data, slong level);
void expr_mul(expr_ptr x, expr_ptr a, expr_ptr b)
{
    _expr_xx_init(x, a, b);
    x->eval = &_expr_mul_eval;
}
void _expr_mul_eval(arb_t z, expr_data_ptr data, slong level)
{
    arb_t u, v;
    slong prec = 1 << level;
    _expr_xx_ptr d = (_expr_xx_ptr) data;
    arb_init(u);
    arb_init(v);
    expr_eval(u, d->a, level);
    expr_eval(v, d->b, level);
    arb_mul(z, u, v, prec);
    arb_clear(u);
    arb_clear(v);
}

void _expr_sub_eval(arb_t z, expr_data_ptr data, slong level);
void expr_sub(expr_ptr x, expr_ptr a, expr_ptr b)
{
    _expr_xx_init(x, a, b);
    x->eval = &_expr_sub_eval;
}
void _expr_sub_eval(arb_t z, expr_data_ptr data, slong level)
{
    arb_t u, v;
    slong prec = 1 << level;
    _expr_xx_ptr d = (_expr_xx_ptr) data;
    arb_init(u);
    arb_init(v);
    expr_eval(u, d->a, level);
    expr_eval(v, d->b, level);
    arb_sub(z, u, v, prec);
    arb_clear(u);
    arb_clear(v);
}

void _expr_div_eval(arb_t z, expr_data_ptr data, slong level);
void expr_div(expr_ptr x, expr_ptr a, expr_ptr b)
{
    _expr_xx_init(x, a, b);
    x->eval = &_expr_div_eval;
}
void _expr_div_eval(arb_t z, expr_data_ptr data, slong level)
{
    arb_t u, v;
    slong prec = 1 << level;
    _expr_xx_ptr d = (_expr_xx_ptr) data;
    arb_init(u);
    arb_init(v);
    expr_eval(u, d->a, level);
    expr_eval(v, d->b, level);
    arb_div(z, u, v, prec);
    arb_clear(u);
    arb_clear(v);
}
