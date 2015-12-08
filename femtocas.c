/*
 * femtocas -- femto-sized computer algebra system
 *
 * Maybe this should be replaced by something like symengine or nemo
 * or at least C++ ...
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
    x->clear = _default_clear;
    x->print = _default_print;
    x->eval = _default_eval;
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
    x->clear = _default_clear;
    x->print = _default_print;
    x->eval = _default_eval;
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



/* make a somewhat complicated expression of three rational numbers */


typedef struct
{
    fmpq a;
    fmpq b;
    fmpq at;
    fmpq bt;
} beta_data_struct;

typedef beta_data_struct * beta_data_ptr;

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

void
beta_expr_init(expr_ptr x, const fmpq_t a, const fmpq_t b, const fmpq_t t)
{
    beta_data_ptr d;

    d = flint_malloc(sizeof(beta_data_struct));
    _fmpq_init_set(&(d->a), a);
    _fmpq_init_set(&(d->b), b);
    fmpq_init(&(d->at));
    fmpq_init(&(d->bt));
    fmpq_mul(&(d->at), a, t);
    fmpq_mul(&(d->bt), b, t);

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
    fmpq_clear(&(d->at));
    fmpq_clear(&(d->bt));
    flint_free(d);
}

void
beta_expr_print(expr_data_ptr data)
{
    beta_data_ptr d = (beta_data_ptr) data;
    flint_printf("beta(");
    flint_printf("a="); fmpq_print(&(d->a));
    flint_printf(" b="); fmpq_print(&(d->b));
    flint_printf(" a*t="); fmpq_print(&(d->at));
    flint_printf(" b*t="); fmpq_print(&(d->bt));
    flint_printf(")");
}

void
beta_expr_eval(arb_t res, expr_data_ptr data, slong level)
{
    slong prec;
    beta_data_ptr d;
    arb_t ax, bx, atx, btx;
    arb_t eatx, ebtx;
    arb_t num, den;

    prec = 1 << level;
    d = (beta_data_ptr) data;

    _arb_init_set_fmpq(ax, &(d->a), prec);
    _arb_init_set_fmpq(bx, &(d->b), prec);
    _arb_init_set_fmpq(atx, &(d->at), prec);
    _arb_init_set_fmpq(btx, &(d->bt), prec);

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

    return 0;
}
