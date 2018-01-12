clear
reset

set title ARG5
set key autotitle columnhead
set offset -0.5

set term pdfcairo noenhanced font "Arial,14"
set output ARG5
#set terminal pdf noenhanced color font 'Helvetica,1'

set style fill solid
set xtics rotate by -25 offset -3,0 nomirror
set xtics font "Arial,11" 
set grid ytics

set ylabel	'Drop Rate in %'
set xlabel	'Backbone Topologies (Rocketfuel)'

print ARG1
print ARG2
print ARG3
print ARG4

plot ARG1  using ($4*100):xtic(1) with histogram title ARG2,\
 ARG3  using ($4*100) with histogram title ARG4,\
 '' using :($4*100+0.6):(sprintf("%d",$2+$3)) with labels font ",10" notitle