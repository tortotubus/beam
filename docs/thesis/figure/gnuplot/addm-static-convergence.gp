if (!exists("csv_file")) csv_file = "addm-static-convergence.csv"
if (!exists("plot_output")) plot_output = "../vector/addm-static-convergence.tex"

set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$\\log{(1/\\Delta s)}$"
set ylabel "$\\log{\\norm{\\vec {\\epsilon}_{\\rm tip}}}$"
set grid
set key off
plot \
  csv_file every ::1 using 6:7 with linespoints title "ADDM"
