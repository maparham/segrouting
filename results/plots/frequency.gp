clear
reset
set term pdfcairo enhanced font "Arial,14"
set output ARG5
set title ARG5

barwidth = 0.3
set boxwidth barwidth absolute
set style fill solid 1.0
set grid ytics

set xlabel "Stack Size"
set ylabel "Frequency"

set xrange [1:4]
set xtics 1,1,4
set format y "%.1l{/Symbol \264}10^%T"
set ytics add ('0' 0)

set offset 0.4,0.4,0,0

print ARG0
print ARG1
print ARG2
print ARG3

plot ARG1 using 4 smooth frequency with boxes title ARG2,\
ARG3 using ($4+barwidth) smooth frequency with boxes title ARG4