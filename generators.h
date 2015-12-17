#ifndef GENERATORS_H
#define GENERATORS_H

#include"expressions.h"


typedef struct generator_reg_node_struct
{
    slong index;
    fmpz * p;
    struct generator_reg_node_struct * next;
} generator_reg_node_struct;
typedef generator_reg_node_struct * generator_reg_node_ptr;
typedef generator_reg_node_struct generator_reg_node_t[1];

typedef struct
{
    generator_reg_node_ptr head;
    generator_reg_node_ptr tail;
    slong size; /* current number of generators */
    slong expressions_len; /* fixed number of expressions */
    int open_for_adding;
} generator_reg_struct;
typedef generator_reg_struct generator_reg_t[1];
typedef generator_reg_struct * generator_reg_ptr;



#ifdef __cplusplus
extern "C" {
#endif


void generator_reg_init(generator_reg_ptr x, slong expressions_len);
slong generator_reg_expressions_len(generator_reg_ptr x);
slong generator_reg_generators_len(generator_reg_ptr x);

void generator_reg_get_matrix(fmpz_mat_t mat, generator_reg_ptr x);
void generator_reg_clear(generator_reg_ptr x);

#ifdef __cplusplus
}
#endif

#endif
