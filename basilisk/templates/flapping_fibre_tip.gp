# Usage: gnuplot -persist -c basilisk/templates/flapping_fibre_tip.gp [path/to/tip.csv]
if (ARGC > 0) { datafile = ARG1 } else { datafile = "huang_output/tip.csv" }
set datafile separator comma
set datafile columnheaders
set xlabel "Time"
set ylabel "Tip position"
set grid
set key outside
# Re-read the CSV every second; press Ctrl+C in the terminal to stop.
while (1) {
    plot datafile using 2:3 with lines title "x tip", \
         datafile using 2:4 with lines title "y tip"
    pause 1
}
