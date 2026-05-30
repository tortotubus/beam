if (!exists("csv_file")) csv_file = "addm-nonlinear-convergence.csv"
if (!exists("plot_output")) plot_output = "../vector/addm-nonlinear-convergence.tex"

set datafile separator comma
set terminal cairolatex pdf size 2.75in,2.25in color colortext font ",7"
set output plot_output

set xlabel "$\\log{(1/\\Delta t)}$" font ",6"
set ylabel "$\\log{\\lVert \\vec x(T) - \\vec x_h(T)\\rVert}$" font ",6"
set grid
set key inside right font ",6" spacing 1.6 width 7 height 1.5

set tics font ",4"
set xlabel offset 0,-0.5
set ylabel offset -1.65,0
set ytics 1.0
set xtics 0.2
set xrange [1.3:2.2]
set yrange [-7.5:-4.5]
plot \
	csv_file every ::1 using ($2 == 1e-07 ? $4 : 1/0):5 with linespoints title "$10^{-7}$", \
	csv_file every ::1 using ($2 == 1e-08 ? $4 : 1/0):5 with linespoints title "$10^{-8}$", \
	csv_file every ::1 using ($2 == 1e-09 ? $4 : 1/0):5 with linespoints title "$10^{-9}$",
