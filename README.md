An experimental dynamic programming implementation
that uses numerical bounds instead of just using double precision
and hoping for the best.

Most input and output uses the json format which is a little nicer
for programmatic access than than using command line flags,
but it is a little trickier to use manually.

Here's a scratch example of using [jq](https://stedolan.github.io/jq/)
to essentially supply a command line parameter (output image filename)
by adding a field to a top level json object on the fly...
`$ jq '.image_filename = "tableau.png"' alignme.json | bin/arbtkf91-image`
![tableau](https://github.com/argriffing/arbtkf91/blob/master/tableau.png)

