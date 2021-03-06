#ifndef TKF91_DP_H
#define TKF91_DP_H

/*
 * This describes the common interface
 * for tkf91 dynamic programming solvers.
 */

#include <stdio.h>
#include <time.h>

#include "flint/flint.h"
#include "flint/fmpz_mat.h"

#include "expressions.h"
#include "dp.h"
#include "tkf91_generator_indices.h"


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Define the data structure that is used for the output.
 * The construction and destruction is defined and performed
 * in the arbtk91.c module; the specific tkf91 dynamic programming
 * solvers themselves do not need to know about those functions.
 */
typedef struct
{
    char * A;
    char * B;
    slong len;
    arb_t log_probability;
    int optimality_flag;
    dp_mat_ptr mat;
} solution_struct;
typedef solution_struct solution_t[1];

void solution_init(solution_t sol, slong aln_maxlen);
void solution_clear(solution_t sol);
void solution_fprint(FILE * file, const solution_t sol);

static __inline__ void
solution_print(const solution_t x)
{
    solution_fprint(stdout, x);
}


/*
 * Define the data structure to be used for miscellaneous request options.
 * The png filename is NULL if no tableau visualization is requested.
 * The trace option is nonzero if the actual alignment is requested (as
 * opposed to requesting just the best score).
 * The rtol option specifies a relative tolerance to be used in the
 * traceback phase of 'float' and 'double' precision dynamic programming;
 * the rtol option is ignored for more sophisticated precision settings.
 */
typedef struct
{
    int trace;
    double rtol;
} request_struct;
typedef request_struct request_t[1];


/* function pointer typedef for the dynamic programming function */
typedef void (*tkf91_dp_fn)(
        solution_t, const request_t,
        fmpz_mat_t, expr_ptr *,
        const tkf91_generator_indices_t,
        const slong *, size_t,
        const slong *, size_t);


#ifdef __cplusplus
}
#endif

#endif
