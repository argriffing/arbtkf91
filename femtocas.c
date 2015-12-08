/*
 * femtocas -- femto-sized computer algebra system
 *
 * Maybe this should be replaced by something like symengine or nemo
 * or actual C++ ...
 *
 * struct
 * x fixed-size array of cached real balls from internal precision doublings
 * x number of such cached precision doublings
 * x pointer to polymorphically defined data
 * x pointer to function to clear the polymorphically defined data
 * x pointer to function that prints given the data
 * x pointer to function that returns real ball given data and internal prec
 *
 * not in struct
 * x details of specialized data structs
 * x other required specialized functions to initialize those structs
 * x generic function dealing with the precision caching
 *
 * NOTE: the interface will use precision levels that are log2 precision bits
 */

#include "flint/flint.h"
#include "flint/fmpq.h"
#include "arb.h"


/* cache at most this many evaluations */
#define EXPR_CACHE_CAP 30

typedef void *expr_data_ptr;

typedef void (*expr_clear_fn)(expr_data_ptr);
typedef void (*expr_print_fn)(expr_data_ptr);
typedef void (*expr_eval_fn)(arb_t, expr_data_ptr, slong);

typedef struct
{
    arb_struct cache[EXPR_CACHE_CAP];
    slong cachesize;
    expr_data_ptr data;
    expr_clear_fn clear;
    expr_print_fn print;
    expr_eval_fn eval;
} expr_struct;

typedef expr_struct expr_t[1];
typedef expr_struct * expr_ptr;


/* helper functions */

void
_fmpq_init_set(fmpq_t z, const fmpq_t x)
{
    fmpq_init(z);
    fmpq_set(z, x);
}

void
_arb_init_set_fmpq(arb_t z, const fmpq_t u, slong prec)
{
    arb_init(z);
    arb_set_fmpq(z, u, prec);
}


/* declare these functions which will be referenced from function pointers */
void beta_expr_clear(expr_data_ptr data);
void beta_expr_print(expr_data_ptr data);
void beta_expr_eval(arb_t res, expr_data_ptr data, slong level);


void
_default_clear(expr_data_ptr p)
{
    flint_printf("clear() not implemented\n");
    abort();
}

void
_default_print(expr_data_ptr p)
{
    flint_printf("print() not implemented\n");
    abort();
}

void
_default_eval(arb_t res, expr_data_ptr p, slong level)
{
    flint_printf("eval() not implemented\n");
    abort();
}

void
expr_init(expr_ptr x)
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

void _expr_i_clear(expr_data_ptr data);
void _expr_i_print(expr_data_ptr data);

void
_expr_i_init(expr_ptr x, const fmpz_t a)
{
    _expr_i_ptr d = flint_malloc(sizeof(_expr_i_struct));
    fmpz_init_set(&(d->a), a);
    expr_init(x);
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

void _expr_q_clear(expr_data_ptr data);
void _expr_q_print(expr_data_ptr data);

void
_expr_q_init(expr_ptr x, const fmpq_t a)
{
    _expr_q_ptr d = flint_malloc(sizeof(_expr_q_struct));
    _fmpq_init_set(&(d->a), a);
    expr_init(x);
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

void _expr_x_clear(expr_data_ptr data);
void _expr_x_print(expr_data_ptr data);

void
_expr_x_init(expr_ptr x, expr_ptr a)
{
    _expr_x_ptr d = flint_malloc(sizeof(_expr_x_struct));
    d->a = a;
    expr_init(x);
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

void _expr_xx_clear(expr_data_ptr data);
void _expr_xx_print(expr_data_ptr data);

void
_expr_xx_init(expr_ptr x, expr_ptr a, expr_ptr b)
{
    _expr_xx_ptr d = flint_malloc(sizeof(_expr_xx_struct));
    d->a = a;
    d->b = b;
    expr_init(x);
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



/* make a somewhat complicated expression of three rational numbers */


typedef struct
{
    fmpq a;
    fmpq b;
    fmpq t;
} beta_data_struct;

typedef beta_data_struct * beta_data_ptr;

void
beta_expr_init(expr_ptr x, const fmpq_t a, const fmpq_t b, const fmpq_t t)
{
    beta_data_ptr d;

    d = flint_malloc(sizeof(beta_data_struct));
    _fmpq_init_set(&(d->a), a);
    _fmpq_init_set(&(d->b), b);
    _fmpq_init_set(&(d->t), t);

    expr_init(x);
    x->data = d;
    x->clear = &beta_expr_clear;
    x->print = &beta_expr_print;
    x->eval = &beta_expr_eval;
}

void
beta_expr_clear(expr_data_ptr data)
{
    beta_data_ptr d = (beta_data_ptr) data;
    fmpq_clear(&(d->a));
    fmpq_clear(&(d->b));
    fmpq_clear(&(d->t));
    flint_free(d);
}

void
beta_expr_print(expr_data_ptr data)
{
    beta_data_ptr d = (beta_data_ptr) data;
    flint_printf("beta(");
    fmpq_print(&(d->a)); flint_printf(", ");
    fmpq_print(&(d->b)); flint_printf(", ");
    fmpq_print(&(d->t));
    flint_printf(")");
}

void
beta_expr_eval(arb_t res, expr_data_ptr data, slong level)
{
    slong prec;
    beta_data_ptr d;
    fmpq_t at, bt;
    arb_t ax, bx, atx, btx;
    arb_t eatx, ebtx;
    arb_t num, den;

    prec = 1 << level;
    d = (beta_data_ptr) data;

    fmpq_init(at);
    fmpq_mul(at, &(d->a), &(d->t));

    fmpq_init(bt);
    fmpq_mul(bt, &(d->b), &(d->t));

    _arb_init_set_fmpq(ax, &(d->a), prec);
    _arb_init_set_fmpq(bx, &(d->b), prec);
    _arb_init_set_fmpq(atx, at, prec);
    _arb_init_set_fmpq(btx, bt, prec);

    arb_init(eatx);
    arb_exp(eatx, atx, prec);

    arb_init(ebtx);
    arb_exp(ebtx, btx, prec);

    arb_init(num);
    arb_sub(num, ebtx, eatx, prec);

    arb_init(den);
    arb_mul(den, bx, ebtx, prec);
    arb_submul(den, ax, eatx, prec);

    arb_div(res, num, den, prec);

    fmpq_clear(at);
    fmpq_clear(bt);

    arb_clear(ax);
    arb_clear(bx);
    arb_clear(atx);
    arb_clear(btx);

    arb_clear(eatx);
    arb_clear(ebtx);

    arb_clear(num);
    arb_clear(den);
}


/* demonstrate usage of the beta expression */


int main()
{
    slong level;
    fmpq_t a, b, t;
    arb_t value;
    expr_t x;

    fmpq_init(a);
    fmpq_init(b);
    fmpq_init(t);

    fmpq_set_si(a, 1, 2);
    fmpq_set_si(b, 3, 4);
    fmpq_set_si(t, 5, 6);

    beta_expr_init(x, a, b, t);
    arb_init(value);

    for (level = 0; level < 8; level++)
    {
        flint_printf("evaluating ");
        expr_print(x);
        flint_printf(" at level %wd:\n", level);
        expr_eval(value, x, level);
        arb_print(value);
        flint_printf("\n\n");
    }

    fmpq_clear(a);
    fmpq_clear(b);
    fmpq_clear(t);
    arb_clear(value);
    expr_clear(x);

    /* do this in a different way */

    fmpq_init(a);
    fmpq_init(b);
    fmpq_init(t);

    fmpq_set_si(a, 1, 2);
    fmpq_set_si(b, 3, 4);
    fmpq_set_si(t, 5, 6);

    fmpq_t at, bt;

    fmpq_init(at);
    fmpq_init(bt);

    fmpq_mul(at, a, t);
    fmpq_mul(bt, b, t);

    expr_t ax, bx;
    expr_t eat, ebt, aeat, bebt;
    expr_t num, den;
    expr_t res;

    expr_fmpq(ax, a);
    expr_fmpq(bx, b);
    expr_exp_fmpq(ebt, bt);
    expr_exp_fmpq(eat, at);
    expr_mul(bebt, bx, ebt);
    expr_mul(aeat, ax, eat);
    expr_sub(num, ebt, eat);
    expr_sub(den, bebt, aeat);
    expr_div(res, num, den);

    for (level = 0; level < 8; level++)
    {
        flint_printf("evaluating ");
        expr_print(res);
        flint_printf(" at level %wd:\n", level);
        expr_eval(value, res, level);
        arb_print(value);
        flint_printf("\n\n");
    }

    expr_clear(ax);
    expr_clear(bx);
    expr_clear(eat);
    expr_clear(ebt);
    expr_clear(aeat);
    expr_clear(bebt);
    expr_clear(num);
    expr_clear(den);
    expr_clear(res);

    return 0;
}
