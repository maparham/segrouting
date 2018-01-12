#! /bin/bash

gnuplot -c failrate.gp \
../data/noflush_TILFA/overall.txt "no flush" \
../data/withflush_TILFA/overall.txt "flush" \
"Packet Drop Rates With TI-LFA, double failures"