#! /bin/bash

TILFA_NOFLUSH="../data/mean_path_length_TILFA_noflush.txt"
TILFA_FLUSH="../data/mean_path_length_TILFA_withflush.txt"
TIMFA_NOFLUSH="../data/mean_path_length_TIMFA_noflush.txt"
TIMFA_FLUSH="../data/mean_path_length_TIMFA_withflush.txt"
rm $TILFA_NOFLUSH $TILFA_FLUSH $TIMFA_NOFLUSH $TIMFA_FLUSH

awk '{sum5[$2] += $3;count4[$2]++;}; END{ for (id in sum5) { print id, sum5[id]/count4[id] } }' < ../data/noflush_TILFA/all.txt >> $TILFA_NOFLUSH
awk '{sum5[$2] += $3;count4[$2]++;}; END{ for (id in sum5) { print id, sum5[id]/count4[id] } }' < ../data/withflush_TILFA/all.txt >> $TILFA_FLUSH
awk '{sum5[$2] += $3;count4[$2]++;}; END{ for (id in sum5) { print id, sum5[id]/count4[id] } }' < ../data/noflush_doubleTILFA/all.txt >> $TIMFA_NOFLUSH
awk '{sum5[$2] += $3;count4[$2]++;}; END{ for (id in sum5) { print id, sum5[id]/count4[id] } }' < ../data/withflush_doubleTILFA/all.txt >> $TIMFA_FLUSH
sleep 1

sort -k 1 $TILFA_NOFLUSH -o $TILFA_NOFLUSH
sort -k 1 $TILFA_FLUSH -o $TILFA_FLUSH
sort -k 1 $TIMFA_NOFLUSH -o $TIMFA_NOFLUSH
sort -k 1 $TIMFA_FLUSH -o $TIMFA_FLUSH
sleep 1
gnuplot -c mean_pathlength.gp \
$TILFA_NOFLUSH "TI-LFA, no flush" \
$TILFA_FLUSH "TI-LFA, flush" \
$TIMFA_NOFLUSH "TI-MFA, no flush" \
$TIMFA_FLUSH "TI-MFA, flush" \
"Average Path Lengths (all 4 algorithms)"
