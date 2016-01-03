#include "jsonutil.h"


static void _fill_sequence_vector(slong *v, const char *str, slong n);

void
_fill_sequence_vector(slong *v, const char *str, slong n)
{
    /* treat N as A following the questionable choice
     * in a reference implementation */
    int i;
    for (i = 0; i < n; i++)
    {
        switch(str[i])
        {
            case 'A' :
            case 'a' :
                v[i] = 0;
                break;
            case 'C' :
            case 'c' :
                v[i] = 1;
                break;
            case 'G' :
            case 'g' :
                v[i] = 2;
                break;
            case 'T' :
            case 't' :
                v[i] = 3;
                break;
            case '-' :
                v[i] = -1;
                break;
            default:
                /* ambiguous nucleotides will be treated as A */
                if (isalpha(str[i]))
                {
                    v[i] = 0;
                }
                else
                {
                    fprintf(stderr, "unrecognized nucleotide ");
                    fprintf(stderr, "ascii %d\n", (int) str[i]);
                    abort();
                }
        }
    }
}


const char *
_json_object_get_string(const json_t *object, const char *key)
{
    json_t *tmp;

    tmp = json_object_get(object, key);
    if (!tmp)
    {
        fprintf(stderr, "error: input does not contain '%s'\n", key);
        abort();
    }
    if (!json_is_string(tmp))
    {
        fprintf(stderr, "error: '%s' is not a json string\n", key);
        abort();
    }
    return json_string_value(tmp);
}


slong *
_json_object_get_sequence(slong *plen, const json_t *object, const char *key)
{
    /* (a, c, g, t, -) -> (0, 1, 2, 3, -1) */
    json_t *tmp;
    const char *value;
    size_t len;
    slong *s;

    tmp = json_object_get(object, key);
    if (!tmp)
    {
        fprintf(stderr, "error: input does not contain '%s'\n", key);
        abort();
    }
    if (!json_is_string(tmp))
    {
        fprintf(stderr, "error: '%s' is not a json string\n", key);
        abort();
    }

    value = json_string_value(tmp);
    len = json_string_length(tmp);

    s = flint_malloc(len * sizeof(slong));
    _fill_sequence_vector(s, value, len);

    *plen = (slong) len;
    return s;
}


slong
_json_object_get_si(const json_t *object, const char *key)
{
    json_t *tmp;
    json_int_t value;

    tmp = json_object_get(object, key);
    if (!tmp)
    {
        fprintf(stderr, "error: input does not contain '%s'\n", key);
        abort();
    }
    if (!json_is_integer(tmp))
    {
        fprintf(stderr, "error: '%s' is not a json integer\n", key);
        abort();
    }

    value = json_integer_value(tmp);
    return (slong) value;
}


void
_json_object_get_fmpq(fmpq_t res, const json_t *object,
        const char *key_n, const char *key_d)
{
    slong n, d;
    n = _json_object_get_si(object, key_n);
    d = _json_object_get_si(object, key_d);
    fmpq_set_si(res, n, d);
}
