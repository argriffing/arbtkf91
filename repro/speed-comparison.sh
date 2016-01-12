SEQUENCES=../../bioinf2015/implementation-team2/doc/sequences

CMD="./TKFLOG_CACHING_ROUND_UP -t"
echo $CMD
pushd ../../bioinf2015/implementation-team2/bin
$CMD
popd
echo

PRECISION=float
CMD="python bench-compare.py --sequences=$SEQUENCES --precision=float"
echo $CMD
$CMD
echo

PRECISION=double
CMD="python bench-compare.py --sequences=$SEQUENCES --precision=double"
echo $CMD
$CMD
echo
