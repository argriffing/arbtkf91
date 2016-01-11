#ifndef RUNJSON_H
#define RUNJSON_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "jansson.h"


/*
 * abuse C as a functional language
 *
 * For the purposes of notation in this file I'm going to pretend that
 * a string homomorphism just means a function that maps a string to a string,
 * and analogously for a json structure homomorphism.
 *
 * I think it is not necessary to have the string homomorphism stuff
 * in this file; it could go into the C file instead.
 */


#ifdef __cplusplus
extern "C" {
#endif

/* a function that maps a string to a string */
typedef char *(*string_hom_fn_t)(void *userdata, const char *data);

/* a function that maps a json struct to a json struct */
typedef json_t *(*json_hom_fn_t)(void *userdata, json_t *data);

/* generic destructor */
typedef void (*clear_fn_t)(void *userdata);

typedef struct
{
    json_hom_fn_t f;
    clear_fn_t clear;
    void *userdata;
} json_hom_struct;
typedef json_hom_struct json_hom_t[1];
typedef json_hom_struct * json_hom_ptr;

typedef struct
{
    string_hom_fn_t f;
    clear_fn_t clear;
    void *userdata;
} string_hom_struct;
typedef string_hom_struct string_hom_t[1];
typedef string_hom_struct * string_hom_ptr;

/* applies a string homomorphism to stdin->stdout */
int run_string_script(string_hom_t hom);

/* creates a string homomorphism induced by a json homomorphism */
string_hom_ptr json_induced_string_hom(json_hom_t hom);

/* induces and applies a string homomorphism to stdin->stdout */
int run_json_script(json_hom_t hom);


#ifdef __cplusplus
}
#endif

#endif
