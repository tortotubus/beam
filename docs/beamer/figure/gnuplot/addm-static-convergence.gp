if (!exists("csv_file")) csv_file = "addm-static-convergence.csv"
if (!exists("plot_output")) plot_output = "../vector/addm-static-convergence.tex"

set datafile separator comma
set terminal cairolatex pdf size 2.75in,2.25in color colortext font ",7"
set output plot_output

set xlabel "$\\log{(1/\\Delta s)}$" font ",6"
set ylabel "$\\log{\\norm{\\vec {\\epsilon}_{\\rm tip}}}$" font ",6"
set grid
set key off
set tics font ",4"

set xlabel offset 0,-0.5
set ylabel offset -1.65,0

set ytics 1.0
set xtics 0.5
plot \
  csv_file every ::1 using 6:7 with linespoints title "ADDM"
