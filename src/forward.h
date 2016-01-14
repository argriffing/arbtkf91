#ifndef FORWARD_H
#define FORWARD_H

#include "dp.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef void *(*forward_strategy_init_t)(void *userdata,
        size_t num);
typedef void (*forward_strategy_clear_t)(void *userdata,
        void *celldata, size_t num);
typedef int (*forward_strategy_visit_t)(void *userdata,
        dp_mat_t mat, slong i, slong j,
        void *curr, void *top, void *diag, void *left);

typedef struct
{
    forward_strategy_init_t init;
    forward_strategy_clear_t clear;
    forward_strategy_visit_t visit;
    size_t sz_celldata;
    void *userdata;
} forward_strategy_struct;
typedef forward_strategy_struct forward_strategy_t[1];
typedef forward_strategy_struct * forward_strategy_ptr;

int dp_forward(dp_mat_t mat, forward_strategy_t strat);


#ifdef __cplusplus
}
#endif

#endif
