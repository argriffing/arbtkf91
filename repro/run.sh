BENCH_DATA="../../../stamatakis/benchMark_data"
TEAM1="../../bioinf2015/implementation-team1/build/mlalign"
TEAM2="../../bioinf2015/implementation-team2/bin/TKFLOG_CACHING_ROUND_UP"
REF="../../../stamatakis/tkf91_scaling/tkf91"

CMD="python repro.py --ref-exe=$REF"
echo $CMD
$CMD

CMD="python repro2.py --exe=$TEAM2 --bench-data=$BENCH_DATA"
echo $CMD
$CMD

CMD="python repro3.py --exe=$TEAM1 --bench-data=$BENCH_DATA"
echo $CMD
$CMD

CMD="python repro4.py --precision=float --bench-data=$BENCH_DATA"
echo $CMD
$CMD

CMD="python bench-analysis.py --precision=float --samples=10 --bench-data=$BENCH_DATA"
echo $CMD
$CMD
