if (!exists("csv_file"))    csv_file = "marchetti-banaei-settling.csv"
if (!exists("plot_output")) plot_output = "../vector/marchetti-banaei-settling.tex"

set datafile separator comma
set terminal cairolatex pdf size 3.2in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$1/(\\gamma r)$"
set ylabel "$\\delta/(L/2)$"

set logscale x
set logscale y

set xrange [10:2000]
set yrange [1e-2:1.2]

set format x "$10^{%T}$"
set format y "$10^{%T}$"

set grid
set size square
set clip two
set key inside bottom right

plot \
  csv_file every ::1 using (strcol(1) eq "bead_spring" ? $2 : 1/0):3 with lines lw 2 dt 1 title "Bead-spring model", \
  csv_file every ::1 using (strcol(1) eq "slender_body" ? $2 : 1/0):3 with lines lw 2 dt 1 title "Slender body", \
  csv_file every ::1 using (strcol(1) eq "present" ? $2 : 1/0):3 with points pt 7 ps 0.9 title "Present results"
