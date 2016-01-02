/*
 * This executable checks the status of a given alignment
 * in the context of a tkf91 model with given parameter values.
 *
 * Possible outcomes:
 * 1) The alignment is shown to be optimal, both in the sense that no other
 *    alignment has greater probability and that it is the canonical
 *    alignment (according to the 'counter-clockwise traceback' criterion)
 *    among all such alignments.
 * 2) The alignment is shown to have the maximum probability
 *    but is not canonical.
 * 3) The alignment is shown to not have the maximum probability.
 * 4) The status could not be determined.
 *
 * This program uses json as the input and output format
 * using the 'jansson' C library to parse the json.
 */

#include "flint/flint.h"

#include "jansson.h"

#include "runjson.h"


/*main(int argc, char *argv[])*/
int
main(void)
{

    flint_cleanup();
    return 0;
}
