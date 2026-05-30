if (!exists("csv_file")) csv_file = "static-sphere-Re-sweep.csv"
if (!exists("plot_output")) plot_output = "../vector/static-sphere-Re-sweep.tex"

set datafile separator comma
set terminal cairolatex pdf size 2.75in,2.25in color colortext font ",7"
set output plot_output

set xlabel "$\\Reynolds_p$"
set ylabel "$C_D$"

set xrange [0:260]
set yrange [0:5.25]

set grid 
set key top right font ",6" spacing 1.6 width 7 height 1.5
set tics font ",6"

set xlabel offset 0,-0.5
set ylabel offset -0.5,0

Cd_SN(Re) = Re > 0 ? 24.0/Re*(1.0 + 0.15*Re**0.687) : 1/0

plot \
	[x=5:275] Cd_SN(x) with lines title "Schiller-Naumann", \
	csv_file every ::1 using 1:2 with points title "DF-IBM"
