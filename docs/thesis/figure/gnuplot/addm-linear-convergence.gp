if (!exists("csv_file")) csv_file = "addm-linear-convergence.csv"
if (!exists("plot_output")) plot_output = "../vector/addm-linear-convergence.tex"

set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$\\log{(1/\\Delta t)}$"
set ylabel "$\\log{\\lVert \\vec x(T) - \\vec x_h(T)\\rVert}$"
set grid
set key off 
plot \
	csv_file every ::1 using 4:5 with linespoints