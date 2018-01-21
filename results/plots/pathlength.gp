clear
reset

set title ARG9
set offset -0.5
set key left Left

#set term pdfcairo noenhanced font "Arial,14"
#set output ARG9
#set terminal pdf noenhanced color font 'Helvetica,1'

set style fill solid
set xtics rotate by -25 offset -2,0 nomirror
#set xtics font "Arial,11" 
set grid ytics

set ylabel	'Path Length'
set xlabel	'Backbone Topologies (Rocketfuel)'

set style data histogram
set style histogram clustered

print ARG1
print ARG3
print ARG5
print ARG7

plot ARG1 using 5:3:xtic(2) with points pt "." title ARG2,\
#	 ARG3  using 2 title ARG4,\
#	 ARG5  using 2 title ARG6,\
#	 ARG7  using 2 title ARG8
