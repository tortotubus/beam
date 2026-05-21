if (!exists("csv_file")) csv_file = "huang-tip-Re-200-KB-10e-4.csv"
if (!exists("plot_output")) plot_output = "../vector/huang-tip-Re-200-KB-10e-4.tex"
if (!exists("gnuplot_source_dir")) gnuplot_source_dir = "."
if (!exists("huang_csv")) huang_csv = sprintf("%s/huang/fig13.csv", gnuplot_source_dir)

set datafile separator comma
set terminal cairolatex pdf size 4.8in,3.2in color colortext font ",10"
set output plot_output

set xlabel "$t$"
set ylabel "$y_{\\rm tip}$"
set xrange [0:50]
set yrange [-0.42:0.42]
set grid
set key top right

set style line 1 lc rgb "#1f77b4" lw 2
set style line 2 lc rgb "#d62728" lw 2 dt 2

plot \
    csv_file every 100::1 using 2:4 with lines ls 1 title "Present", \
    huang_csv every ::1 using 1:2 with lines ls 2 title "Huang et al."
