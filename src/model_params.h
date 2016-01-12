#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include "flint/flint.h"
#include "flint/fmpq.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
    fmpq_t lambda;
    fmpq_t mu;
    fmpq_t tau;
    fmpq pi[4];
} model_params_struct;
typedef model_params_struct model_params_t[1];
typedef model_params_struct * model_params_ptr;

void model_params_init(model_params_t p);
void model_params_clear(model_params_t p);
void model_params_print(const model_params_t p);

/* complain on stderr and return a nonzero value if parameters are invalid */
int model_params_validate(const model_params_t p);


#ifdef __cplusplus
}
#endif

#endif
