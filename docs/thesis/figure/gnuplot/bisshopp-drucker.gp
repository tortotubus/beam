if (!exists("csv_file")) csv_file = "bisshopp-drucker.csv"
if (!exists("plot_output")) plot_output = "../vector/bisshopp-drucker.tex"

set datafile separator comma
set terminal cairolatex pdf size 3.2in,3.2in color colortext font ",10"
set output plot_output

set xlabel ""
set ylabel "$PL^2/B$"
set xtics 0.0,0.1,1.00
set ytics 0.0,1.0,10.0
set xrange [0:1]
set yrange [0:10]
set grid
set size square
set key inside left 
plot \
  csv_file every ::1 using 1:3 with lines title "$x_{\\rm tip}/ L$" smooth csplines, \
  csv_file every ::1 using 2:3 with lines title "$y_{\\rm tip}/ L$" smooth csplines
