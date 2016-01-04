#include <stdio.h>

#include "jansson.h"

#include "runjson.h"


char *jsonwrap(void *userdata, const char *s_in);

char *jsonwrap(void *userdata, const char *s_in)
{
    json_hom_ptr p = userdata;

    json_error_t error;
    const size_t flags = 0;
    json_t *j_in;
    json_t *j_out;
    char *s_out;
    char name[] = "json wrapper";

    /* convert input string to json using jansson */
    j_in = json_loads(s_in, flags, &error);
    if (!j_in) {
        fprintf(stderr, "json error in %s: on json line %d: %s\n",
            name, error.line, error.text);
        return NULL;
    }

    /* call the underlying function to get a new jansson object */
    j_out = p->f(p->userdata, j_in);
    if (!j_out) {
        /*
        fprintf(stderr, "error: failed to get json object output\n");
        return NULL;
        */
    }

    /* convert the output jansson object to a string */
    if (j_out)
    {
        s_out = json_dumps(j_out, flags);
        if (!s_out) {
            fprintf(stderr, "error: failed to dump ");
            fprintf(stderr, "the json object to a string\n");
            return NULL;
        }
    }
    else
    {
        s_out = NULL;
    }

    /* attempt to free the input and output json objects */
    json_decref(j_in);
    json_decref(j_out);

    /* return the output string which must be freed by the caller */
    return s_out;
}

string_hom_ptr
json_induced_string_hom(json_hom_t hom)
{
    string_hom_ptr p = malloc(sizeof(string_hom_struct));
    p->userdata = hom;
    p->clear = NULL;
    p->f = jsonwrap;
    return p;
}




char *fgets_dynamic(FILE *stream);

/*
 * Reads unlimited input into a string that must be freed by the caller.
 * This is just a utility function.
 */
char *fgets_dynamic(FILE *stream)
{
  int size = 0;
  int capacity = 20;
  char *s = calloc(capacity, sizeof(*s));
  char *tail = s;

  /* keep reading from the command line until EOF */
  while (fgets(tail, capacity - size, stream)) {
    /* fprintf(stderr, "debug: {size: %d, capacity: %d}\n", size, capacity);
     */
    size += strlen(tail);
    if (size == capacity-1) {
      capacity <<= 1;
      s = realloc(s, capacity * sizeof(*s));
      if (!s) {
        fprintf(
            stderr,
            "error: fgets_dynamic: failed to reallocate %d\n",
            capacity);
        return 0;
      }
    }
    tail = s + size;
  }

  /* return the newly allocated and filled string */
  return s;
}



int
run_string_script(string_hom_t hom)
{
  char *s_in = 0;
  char *s_out = 0;
  char name[] = "string->string wrapper";

  /* read the string from stdin until we find eof */
  s_in = fgets_dynamic(stdin);
  if (!s_in) {
    fprintf(stderr, "error: %s: failed to read string from stdin\n", name);
    return -1;
  }

  /* call the string interface wrapper */
  s_out = hom->f(hom->userdata, s_in);
  /*
  if (!s_out) {
    fprintf(stderr, "error: %s: failed to get a response\n", name);
    return -1;
  }
  */

  /* free the input string */
  free(s_in);

  /* write the string to stdout */
  /* free the string allocated by the script function */
  if (s_out)
  {
      puts(s_out);
      free(s_out);
  }

  /* return zero if no error */
  return 0;
}


int
run_json_script(json_hom_t hom)
{
    string_hom_ptr p = json_induced_string_hom(hom);
    int result = run_string_script(p);
    free(p);
    return result;
}
