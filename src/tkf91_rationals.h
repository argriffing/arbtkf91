#ifndef TKF91_RATIONALS_H
#define TKF91_RATIONALS_H

#include "flint/flint.h"
#include "flint/fmpq.h"


typedef struct
{
    fmpq_t lambda;
    fmpq_t mu;
    fmpq_t tau;
    fmpq pi[4];

    fmpq qi[4];
    fmpq_t negdt;
    fmpq_t lambda_div_mu;
    fmpq_t one_minus_lambda_div_mu;
    fmpq_t beta_exponent; /* (lambda - mu) * tau */
    fmpq_t neg_mu_tau;

} tkf91_rationals_struct;
typedef tkf91_rationals_struct tkf91_rationals_t[1];
typedef tkf91_rationals_struct * tkf91_rationals_ptr;


#ifdef __cplusplus
extern "C" {
#endif


void tkf91_rationals_init(tkf91_rationals_t r,
        const fmpq_t lambda, const fmpq_t mu, const fmpq_t tau,
        const fmpq * pi);
void tkf91_rationals_clear(tkf91_rationals_t r);


#ifdef __cplusplus
}
#endif

#endif
