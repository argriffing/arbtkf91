#ifndef VIS_H
#define VIS_H

#include "breadcrumbs.h"


#ifdef __cplusplus
extern "C" {
#endif


int write_tableau_image(
        char * filename, const breadcrumb_mat_t mat, char * title);


#ifdef __cplusplus
}
#endif

#endif
