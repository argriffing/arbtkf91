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


json_t *run(void * userdata, json_t *j_in);

json_t *run(void * userdata, json_t *j_in)
{
    if (userdata)
    {
        fprintf(stderr, "error: unexpected userdata\n");
        abort();
    }

    json_t *root;
    json_t *j_out;

    root = j_in;

    if (!json_is_array(root))
    {
        fprintf(stdout, "root is not an array\n");
    }
    else
    {
        fprintf(stdout, "root is an array!\n");
    }

    /* j_out = json_pack("i", 42); */
    j_out = json_pack("[ssb]", "foo", "bar", 1);
    return j_out;
}


int main(void)
{
    json_hom_t hom;
    hom->userdata = NULL;
    hom->clear = NULL;
    hom->f = run;
    int result = run_json_script(hom);

    flint_cleanup();
    return result;
}
