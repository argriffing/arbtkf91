#ifndef UNUSED_H
#define UNUSED_H

/*
 * We want to maintain a consistent interface to allow using function pointers,
 * but sometimes parameters will be unused. According to the internet,
 * this #define is a relatively good way to mark the unused parameters.
 */
#define UNUSED(x) (void)(x)

#endif
