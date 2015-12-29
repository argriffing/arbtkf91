#ifndef PRINTUTIL_H
#define PRINTUTIL_H

#include <stdio.h>
#include <time.h>


#ifdef __cplusplus
extern "C" {
#endif

static __inline__ void
_print_elapsed_time(clock_t diff)
{
    int usec = (diff * 1000 * 1000) / CLOCKS_PER_SEC;
    int msec = usec / 1000;
    int sec = msec / 1000;
    printf("time taken %d seconds %d milliseconds %d microseconds.\n",
            sec, msec%1000, usec%1000);
}

#ifdef __cplusplus
}
#endif

#endif
