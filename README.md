An experimental tkf91 dynamic programming implementation
that uses numerical bounds instead of just using double precision
and hoping for the best.

The four scripts arbtkf91-{align, bench, check, image}
align, benchmark, check alignment optimality,
and generate tableau visualizations respectively.

Most input and output uses the json format which is nicer
for programmatic access than using command line flags,
but it is a little trickier to use manually.

Examples
--------

Some of these examples use [jq](https://stedolan.github.io/jq/)
to filter the json inputs or outputs.

`$ cat in.json`

```javascript
{
    "precision" : "float",
    "rtol" : 0.0,
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

`$ arbtkf91-align < in.json | jq '. | {a: .sequence_a, b: .sequence_b}'`

```javascript
{
  "a": "ACGACTAGTCA-GC-TACG-AT-CGA-CT-C-ATTCAACTGACTGACA-TCGACTTA",
  "b": "A-GAG-AGTAATGCATACGCATGC-ATCTGCTATTC---TG-CTG-CAGTGG--T-A"
}
```

`$ jq '.image_mode="full" | .image_filename="tableau.png"' in.json | bin/arbtkf91-image`

![tableau](https://github.com/argriffing/arbtkf91/blob/master/tableau.png)


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

`$ arbtkf91-align < in.json | arbtkf91-check`

```javascript
{"alignment_is_optimal": "yes", "alignment_is_canonical": "yes", "number_of_optimal_alignments": "56"}
```

These examples show that float precision is not enough
for reasonable sized alignments.

`$ jq '.precision="float" | .rtol=3e-7' fails-high-tolerance.json | bin/arbtkf91-align | bin/arbtkf91-check`
```json
{"alignment_is_optimal": "no", "number_of_optimal_alignments": "32332559983411306514373848819744479641600", "alignment_is_canonical": "no"}
```
`$ jq '.precision="float" | .rtol=3e-7' needs-high-tolerance.json | bin/arbtkf91-align | bin/arbtkf91-check`
```json
{"alignment_is_optimal": "yes", "alignment_is_canonical": "no", "number_of_optimal_alignments": "9442009665687106671596887819668655696812107909520913524435008004699019468288819200000000000000"}
```



`$ jq '.image_mode="simple" | .image_filename="needs-high.tableau.png"' needs-high-tolerance.json | bin/arbtkf91-image`

![tableau](https://github.com/argriffing/arbtkf91/blob/master/needs-high.tableau.png)
