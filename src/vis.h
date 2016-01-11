#ifndef VIS_H
#define VIS_H

#include "breadcrumbs.h"


#ifdef __cplusplus
extern "C" {
#endif


int write_tableau_image(const char * filename,
        const breadcrumb_mat_t mat, const char * title);

int write_simple_tableau_image(const char * filename,
        const breadcrumb_mat_t mat, const char * title);


#ifdef __cplusplus
}
#endif

#endif
