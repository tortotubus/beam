if (!exists("csv_file")) csv_file = "addm-nonlinear-convergence.csv"
if (!exists("plot_output")) plot_output = "../vector/addm-nonlinear-convergence.tex"

set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$\\log{(1/\\Delta t)}$"
set ylabel "$\\log{\\lVert \\vec x(T) - \\vec x_h(T)\\rVert}$"
set grid
set key inside right 
plot \
	csv_file every ::1 using ($2 == 1e-07 ? $4 : 1/0):5 with linespoints title "$\\epsilon = 1e-07$", \
	csv_file every ::1 using ($2 == 1e-08 ? $4 : 1/0):5 with linespoints title "$\\epsilon = 1e-08$", \
	csv_file every ::1 using ($2 == 1e-09 ? $4 : 1/0):5 with linespoints title "$\\epsilon = 1e-09$",