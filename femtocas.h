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
 * NOTE: the interface uses precision levels that are log2 precision bits
 */

#ifndef FEMTOCAS_H
#define FEMTOCAS_H

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"

/*
 * We want to maintain a consistent interface to allow using function pointers,
 * but sometimes parameters will be unused. According to the internet,
 * this #define is a relatively good way to mark the unused parameters.
 */
#define UNUSED(x) (void)(x)

/* cache at most this many evaluations */
#define EXPR_CACHE_CAP 30


/*
 * Forward declaration of expr_struct_tag to allow declaration of expr_ptr.
 * Technically this section could be combined with the next section,
 * but keeping them separate maintains some design flexibility.
 * This section will always stay in the header file, whereas the section
 * that defines the expr_struct could possibly be moved into the C file.
 */
struct expr_struct_tag;
typedef struct expr_struct_tag * expr_ptr;

/*
 * We will define the expr_struct in the header file as well.
 * The idea is to let users create these on the stack,
 * rather than requiring that they be created on the heap.
 */
typedef void *expr_data_ptr;
typedef void (*expr_clear_fn)(expr_data_ptr);
typedef void (*expr_print_fn)(expr_data_ptr);
typedef void (*expr_eval_fn)(arb_t, expr_data_ptr, slong);
typedef struct expr_struct_tag
{
    arb_struct cache[EXPR_CACHE_CAP];
    slong cachesize;
    expr_data_ptr data;
    expr_clear_fn clear;
    expr_print_fn print;
    expr_eval_fn eval;
    void *userdata;
} expr_struct;
typedef expr_struct expr_t[1];


#ifdef __cplusplus
extern "C" {
#endif


/* create the expression */

/* unary(fmpz) */
void expr_fmpz(expr_ptr x, const fmpz_t a);
void expr_log_fmpz(expr_ptr x, const fmpz_t a);

/* unary(fmpq) */
void expr_fmpq(expr_ptr x, const fmpq_t a);
void expr_log_fmpq(expr_ptr x, const fmpq_t a);
void expr_exp_fmpq(expr_ptr x, const fmpq_t a);
void expr_complement_fmpq(expr_ptr x, const fmpq_t a);
void expr_neg_fmpq(expr_ptr x, const fmpq_t a);

/* unary(expr) */
void expr_exp(expr_ptr x, expr_ptr a);
void expr_neg(expr_ptr x, expr_ptr a);
void expr_log(expr_ptr x, expr_ptr a);
void expr_log1p(expr_ptr x, expr_ptr a);
void expr_log1m(expr_ptr x, expr_ptr a);
void expr_complement(expr_ptr x, expr_ptr a);

/* binary(expr, expr) */
void expr_add(expr_ptr x, expr_ptr a, expr_ptr b);
void expr_mul(expr_ptr x, expr_ptr a, expr_ptr b);
void expr_sub(expr_ptr x, expr_ptr a, expr_ptr b);
void expr_div(expr_ptr x, expr_ptr a, expr_ptr b);


/* do things with the created expression */

void expr_print(expr_ptr x);
void expr_eval(arb_t res, expr_ptr x, slong level);


/* finally, clear the memory of the expression */

void expr_clear(expr_ptr x);


#ifdef __cplusplus
}
#endif

#endif
