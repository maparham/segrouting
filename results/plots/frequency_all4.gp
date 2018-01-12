clear
reset
set term pdfcairo enhanced font "Arial,14"
set output ARG9
set title ARG9

barwidth = 0.13
set boxwidth barwidth absolute
set style fill solid 1.0 noborder
set grid ytics
set key Left

set xlabel "Stack Size"
set ylabel "Frequency"

set xrange [1:3]
set xtics 1,1,3 nomirror
set format y "%.1l{/Symbol \264}10^%T"
set ytics add ('0' 0)

set offset 0.4,0.4,0,0

print ARG0
print ARG1
print ARG2
print ARG3

plot ARG1 using 4 smooth frequency with boxes title ARG2,\
ARG3 using ($4+barwidth) smooth frequency with boxes title ARG4,\
ARG5 using ($4+2*barwidth) smooth frequency with boxes title ARG6,\
ARG7 using ($4+3*barwidth) smooth frequency with boxes title ARG8