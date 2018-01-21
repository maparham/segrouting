clear
reset
set term pdfcairo enhanced font "Arial,14"
set output ARG9
set title ARG9

EVERY=1
barwidth = 0.1
set boxwidth barwidth absolute
set style fill solid 1.0 noborder
set grid ytics
set key Left

set xlabel "Stack Size"
set ylabel "Frequency in %"

set xrange [0.85:3.36]
set xtics 1,1,3 nomirror left offset 4.5,0
#set format y "%.1l{/Symbol \264}10^%T"
set ytics add ('0' 0)

#set offset 0.2,0.4,0,0

print ARG0
print ARG1
print ARG2
print ARG3



stats ARG1 using 4 name "A"
stats ARG3 using 4 name "B"
stats ARG5 using 4 name "C"
stats ARG7 using 4 name "D"

plot ARG1 every EVERY using 4:(100.0/A_records) smooth frequency with boxes title ARG2,\
ARG3 every EVERY using ($4+barwidth):(100.0/B_records) smooth frequency with boxes title ARG4,\
ARG5 every EVERY using ($4+2*barwidth+0.05):(100.0/C_records) smooth frequency with boxes title ARG6,\
ARG7 every EVERY using ($4+3*barwidth+0.05):(100.0/D_records)smooth frequency with boxes title ARG8