#ifndef JSONUTIL_H
#define JSONUTIL_H


#include "jansson.h"

#include "flint/flint.h"
#include "flint/fmpq.h"



#ifdef __cplusplus
extern "C" {
#endif

const char * _json_object_get_string(const json_t *object, const char *key);

slong *
_json_object_get_sequence(slong *plen, const json_t *object, const char *key);

slong _json_object_get_si(const json_t *object, const char *key);

void _json_object_get_fmpq(fmpq_t res, const json_t *object,
        const char *key_n, const char *key_d);



#ifdef __cplusplus
}
#endif

#endif
