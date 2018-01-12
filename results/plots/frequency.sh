#! /bin/bash

gnuplot -c frequency_all4.gp \
../data/noflush_TILFA/all.txt "TI-LFA, no flush" \
../data/withflush_TILFA/all.txt  "TI-LFA, flush" \
../data/noflush_doubleTILFA/all.txt "TI-MFA, no flush" \
../data/withflush_doubleTILFA/all.txt "TI-MFA, flush" \
"Frequency of Stack Sizes (all 4 algorithms)"