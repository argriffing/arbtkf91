Introduction
------------

The arbtkf91 program tries to compute the maximum probability tkf91 alignment
for a given sequence pair and set of parameter values.
The numerical error is controlled in such a way that for many inputs
it is possible to verify that the alignment is optimal.

The command line program uses json on stdin and stdout as follows.

`$ arbtkf91-align < examples/in.json`

```javascript
{
    "parameters" :
    {
        "pa" : {"num" : 25, "denom" : 100},
        "pc" : {"num" : 25, "denom" : 100},
        "pg" : {"num" : 25, "denom" : 100},
        "pt" : {"num" : 25, "denom" : 100},
        "lambda" : {"num" : 1, "denom" : 1},
        "mu" : {"num" : 2, "denom" : 1},
        "tau" : {"num" : 1, "denom" : 10}
    },
    "sequence_a": "ACGACTAGTCA-GC-TACG-AT-CGA-CT-C-ATTCAACTGACTGACA-TCGACTTA",
    "sequence_b": "A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A"
}
```

`$ cat examples/in.json`

```javascript
{
    "precision" : "mag",
    "rtol" : 0,
    "parameters" :
    {
        "pa" : {"num" : 25, "denom" : 100},
        "pc" : {"num" : 25, "denom" : 100},
        "pg" : {"num" : 25, "denom" : 100},
        "pt" : {"num" : 25, "denom" : 100},
        "lambda" : {"num" : 1, "denom" : 1},
        "mu" : {"num" : 2, "denom" : 1},
        "tau" : {"num" : 1, "denom" : 10}
    },
    "sequence_a" : "ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA",
    "sequence_b" : "AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA"
}
```


Requirements
------------

The arbtkf91 code has been tested only on Linux.

It depends on these C libraries:
 * [arb](https://github.com/fredrik-johansson/arb)
   -- C library for arbitrary-precision interval arithmetic
 * [flint2](https://github.com/wbhart/flint2)
   -- <b>F</b>ast <b>Li</b>brary for <b>N</b>umber <b>T</b>heory
 * [gmp](https://gmplib.org/)
   -- The GNU Multiple Precision Arithmetic Library
 * [libpng](http://www.libpng.org/pub/png/libpng.html)
   -- The official PNG reference library
 * [jansson](https://github.com/akheron/jansson)
   -- C library for encoding, decoding and manipulating JSON data

The tests depend on a few Python packages:
 * [numpy](https://github.com/numpy/numpy)
   -- A package for scientific computing with Python
 * [biopython](https://github.com/biopython/biopython)
   -- Python tools for computational molecular biology
 * [nose](https://nose.readthedocs.org)
   -- A package for unit testing

Scripts may use [jq](https://stedolan.github.io/jq/) for json filtering.


Installation
------------

Something like the usual autotools installation commands should
work if you are lucky:

```shell
$ ./configure CPPFLAGS='-I/path/to/include/flint'
$ make
$ make check
$ make install
```

The extra CPPFLAGS path is due to the
[idiosyncratic](https://github.com/fredrik-johansson/arb/issues/24)
way that arb includes the flint2 headers.

To use a configuration tuned to your specific machine architecture,
use a `configure` command like the following:

```shell
$ ./configure CFLAGS='-O3 -march=native -ffast-math' CPPFLAGS='-I/path/to/include/flint'
```


Testing
-------

Some tests are run during the `$ make check` step of the installation.
After installation, additional python test scripts can be run
with the `$ nosetests -v` command.


Examples
--------

These example commands assume that the current directory is `examples`.

### align

`$ arbtkf91-align < in.json | jq '. | {a: .sequence_a, b: .sequence_b}'`


```javascript
{
  "a": "ACGACTAGTCA-GC-TACG-AT-CGA-CT-C-ATTCAACTGACTGACA-TCGACTTA",
  "b": "A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A"
}
```


### bench

`$ jq '.samples=10 | .precision="float"' in1k.json | arbtkf91-bench | jq '. | .elapsed_ticks'`
```javascript
[5206, 4800, 3439, 3425, 3425, 3412, 3434, 3420, 3413, 3410]
```

`$ jq '.samples=10 | .precision="double"' in1k.json | arbtkf91-bench | jq '. | .elapsed_ticks'`
```javascript
[7957, 7407, 4340, 4354, 4357, 4315, 4350, 4323, 4336, 4324]
```

### check

`$ arbtkf91-align < in.json | arbtkf91-check`
```javascript
{
    "alignment_is_optimal": true,
    "alignment_is_canonical": true
}
```


`$ jq '.precision="float" | .rtol=3e-7' fails-high-tolerance.json | arbtkf91-align | arbtkf91-check`
```javascript
{
    "alignment_is_optimal": false,
    "alignment_is_canonical": false
}
```

`$ jq '.precision="float" | .rtol=3e-7' needs-high-tolerance.json | arbtkf91-align | arbtkf91-check`
```javascript
{
    "alignment_is_optimal": true,
    "alignment_is_canonical": false
}
```


### count

`$ arbtkf91-count < in.json`
```javascript
{"number_of_optimal_alignments": "56"}
```


### image

`$ jq '.image_mode="full" | .image_filename="tableau.png"' in.json | arbtkf91-image`

![tableau](https://github.com/argriffing/arbtkf91/blob/master/examples/tableau.png)


`$ jq '.image_mode="simple" | .image_filename="needs-high.tableau.png"' needs-high-tolerance.json | arbtkf91-image`

![tableau](https://github.com/argriffing/arbtkf91/blob/master/examples/needs-high.tableau.png)
