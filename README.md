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
    "sequence_b": "A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A",
    "verified": true
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

For alignments and parameter settings that are trickier than the
reference data settings used for benchmarking,
the computed bounds are not always tight enough.
The example below creates a bunch of debugging spam
because the implementation does not yet adjust the precision on the fly;
I have only seen this situation occur for artificially constructed
inputs, so I have not yet bothered to implement an adaptive search.

`$ arbtkf91-align < difficult.json`

```
expected two vectors to be equal, but they aren't
11  18 7 19 2 6 2 5 -28 5 8 8 : (-71250854651360737784577975918036393052090092135798598252541061727090738870129 * 2^-249) +/- (810343202 * 2^-273)
11  18 7 19 2 6 3 4 -29 5 7 8 : (-71250854564012816899681735519152050743298744969918459614906269546550436053621 * 2^-249) +/- (839718047 * 2^-273)
overlap? 0
fail: want3 m2 vs m0
i=21 j=22
{"parameters": {"pa": {"num": 10000001, "denom": 40000010}, "tau": {"num": 1, "denom": 10}, "pc": {"num": 10000002, "denom": 40000010}, "pt": {"num": 10000004, "denom": 40000010}, "pg": {"num": 10000003, "denom": 40000010}, "mu": {"num": 2, "denom": 1}, "lambda": {"num": 1, "denom": 1}}, "sequence_a": "ACGACTAGTCA-GC-TACG-AT-CGA-CT-C-ATTCAACTGACTGACA-TCGACTTA", "sequence_b": "A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A", "verified": false}
```

Sometimes this can be worked around by using a higher, slower,
hardcoded precision setting.
In this case the `arb256` setting
verifies that the solution is the canonical optimal solution.

`$ jq '. | .precision="arb256"' difficult.json | arbtkf91-align`

```
{"sequence_b": "A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A", "parameters": {"pa": {"num": 10000001, "denom": 40000010}, "pc": {"num": 10000002, "denom": 40000010}, "pg": {"num": 10000003, "denom": 40000010}, "pt": {"num": 10000004, "denom": 40000010}, "tau": {"num": 1, "denom": 10}, "lambda": {"num": 1, "denom": 1}, "mu": {"num": 2, "denom": 1}}, "verified": true, "sequence_a": "ACGACTAGTCA-GC-TACG-AT-CGA-CT-C-ATTCAACTGACTGACA-TCGACTTA"}
```



### bench

`$ jq '.samples=10 | .precision="arb256"' in.json | arbtkf91-bench | jq '. | .elapsed_ticks'`

```javascript
[
  16288,
  15987,
  16082,
  16036,
  15957,
  15995,
  15999,
  15949,
  16047,
  16047
]
```

### check

`$ arbtkf91-align < in.json | jq '. | del(.verified)' | arbtkf91-check`
```javascript
{
    "alignment_is_optimal": "yes",
    "alignment_is_canonical": "yes",
    "number_of_optimal_alignments": "56"
}
```

`$ jq '.precision="float" | .rtol=3e-7' fails-high-tolerance.json | arbtkf91-align | jq '. | del(.verified)' | arbtkf91-check`
```json
{
    "alignment_is_optimal": "no",
    "alignment_is_canonical": "no",
    "number_of_optimal_alignments": "32332559983411306514373848819744479641600"
}
```

`$ jq '.precision="float" | .rtol=3e-7' needs-high-tolerance.json | arbtkf91-align | jq '. | del(.verified)' | arbtkf91-check`
```json
{
    "alignment_is_optimal": "yes",
    "alignment_is_canonical": "no",
    "number_of_optimal_alignments": "9442009665687106671596887819668655696812107909520913524435008004699019468288819200000000000000"
}
```


### image

`$ jq '.image_mode="full" | .image_filename="tableau.png"' in.json | arbtkf91-image`

![tableau](https://github.com/argriffing/arbtkf91/blob/master/examples/tableau.png)


`$ jq '.image_mode="simple" | .image_filename="needs-high.tableau.png"' needs-high-tolerance.json | arbtkf91-image`

![tableau](https://github.com/argriffing/arbtkf91/blob/master/examples/needs-high.tableau.png)
