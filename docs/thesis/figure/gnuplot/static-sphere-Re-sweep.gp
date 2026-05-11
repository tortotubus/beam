if (!exists("csv_file")) csv_file = "static-sphere-Re-sweep.csv"
if (!exists("plot_output")) plot_output = "../vector/static-sphere-Re-sweep.tex"

set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$\\Reynolds_p$"
set ylabel "$C_D$"

set xrange [0:260]
set yrange [0:5.25]

set grid 

Cd_SN(Re) = Re > 0 ? 24.0/Re*(1.0 + 0.15*Re**0.687) : 1/0

plot \
	[x=5:275] Cd_SN(x) with lines title "Schiller-Naumann", \
	csv_file every ::1 using 1:2 with points title "DF-IBM $d_p/\\Delta x = 25.6$"
