python repro.py --ref-exe=../../../stamatakis/tkf91_scaling/tkf91
first example...
/path/to/exe --pa 0.27 --pc 0.24 --pg 0.26 --pt 0.23 --lambda 1 --mu 2 --tau 0.1 --sequence-1 ACGACTAGTCAGCTACGATCGACTCATTCAACTGACTGACATCGACTTA --sequence-2 AGAGAGTAATGCATACGCATGCATCTGCTATTCTGCTGCAGTGGTA
It seemed to work.

second example...
/path/to/exe --pa 0.27 --pc 0.24 --pg 0.26 --pt 0.23 --lambda 1 --mu 2 --tau 0.1 --sequence-1 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --sequence-2 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
It seemed to work.

third example...
/path/to/exe --pa 0.27 --pc 0.24 --pg 0.26 --pt 0.23 --lambda 1 --mu 2 --tau 0.1 --sequence-1 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --sequence-2 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
It seemed to fail with the following error message:
tkf91: dp.c:604: tkf91_dp: Assertion `dir >= 0' failed.

python repro2.py --exe=../../bioinf2015/implementation-team2/bin/TKFLOG_CACHING_ROUND_UP --bench-data=../../../stamatakis/benchMark_data
RBP3_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 12

vWF_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 14

RAG2_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 18

cytb_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 13

RAG1_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 12

BDNF_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 13

python repro3.py --exe=../../bioinf2015/implementation-team1/build/mlalign --bench-data=../../../stamatakis/benchMark_data
RBP3_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 11

vWF_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 16

RAG2_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 18

cytb_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 12

RAG1_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 10

BDNF_unaligned.fas
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 13

python repro4.py --precision=float --bench-data=../../../stamatakis/benchMark_data
RBP3_unaligned.fas
precision=float rtol=0.0
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 5

vWF_unaligned.fas
precision=float rtol=0.0
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 3

RAG2_unaligned.fas
precision=float rtol=0.0
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 16

cytb_unaligned.fas
precision=float rtol=0.0
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 2

RAG1_unaligned.fas
precision=float rtol=0.0
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 12

BDNF_unaligned.fas
precision=float rtol=0.0
total number of alignments: 45
number of optimal alignments: 45
number that are also canonical: 16

