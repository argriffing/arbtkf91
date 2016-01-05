#ifndef JSON_MODEL_PARAMS_H
#define JSON_MODEL_PARAMS_H

#include "flint/flint.h"
#include "flint/fmpq.h"

#include "jansson.h"

#include "model_params.h"


#ifdef __cplusplus
extern "C" {
#endif


int _json_get_fmpq(fmpq_t x, json_t * root);
int _json_get_fmpq_ex(fmpq_t x, json_t * root,
        json_error_t *perror, size_t flags);

int _json_get_model_params(model_params_t p, json_t * root);
int _json_get_model_params_ex(model_params_t p, json_t * root,
        json_error_t *perror, size_t flags);


#ifdef __cplusplus
}
#endif

#endif
