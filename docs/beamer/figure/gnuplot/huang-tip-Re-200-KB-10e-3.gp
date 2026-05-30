if (!exists("csv_file")) csv_file = "huang-tip-Re-200-KB-10e-4.csv"
if (!exists("plot_output")) plot_output = "../vector/huang-tip-Re-200-KB-10e-4.tex"
if (!exists("gnuplot_source_dir")) gnuplot_source_dir = "."
if (!exists("huang_csv")) huang_csv = sprintf("%s/huang/fig13.csv", gnuplot_source_dir)

set datafile separator comma
set terminal cairolatex pdf size 5in,1.7in color colortext font ",7"
set output plot_output

set xlabel "$t$"
set ylabel "$y_{\\rm tip}$"
set xrange [0:50]
set yrange [-0.4:0.4]
set ytics 0.2
set grid
set key top right font ",6" spacing 1.2 width 7 height 1.5 box opaque fc rgb "white"
set tics font ",6"
set xlabel offset 0,0.3
set ylabel offset 1.2,0

set style line 1 lc rgb "#1f77b4" lw 5
set style line 2 lc rgb "#d62728" lw 5 dt 2

plot \
    csv_file every 100::1 using 2:4 with lines ls 1 title "Present", \
    huang_csv every ::1 using 1:2 with lines ls 2 title "Huang et al."
