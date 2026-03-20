# Usage:
#   gnuplot -e "file='forces.txt'" plot_forces.gp
# or edit the filename below directly.

if (!exists("file")) file = "/home/colive/Git/gitlab-math/colive/beam/build/debug/basilisk/examples/serial/cylinderstats.txt"

set datafile separator whitespace
set key outside
set grid
set xlabel "t"

set terminal qt 0

# Plot all quantities on one figure
set title "Force and coefficient history"
plot \
  file using 1:2 with lines linewidth 2 title "Fx", \
  file using 1:3 with lines linewidth 2 title "Fy", \
  file using 1:4 with lines linewidth 2 title "Cd", \
  file using 1:5 with lines linewidth 2 title "Cf"

pause -1 "Press Enter to continue"

# Separate plots
set multiplot layout 2,2 rowsfirst title "Force and coefficient history"

set title "Fx vs t"
set xlabel "t"
set ylabel "Fx"
plot file using 1:2 with lines linewidth 2 title "Fx"

set title "Fy vs t"
set xlabel "t"
set ylabel "Fy"
plot file using 1:3 with lines linewidth 2 title "Fy"

set title "Cd vs t"
set xlabel "t"
set ylabel "Cd"
plot file using 1:4 with lines linewidth 2 title "Cd"

set title "Cf vs t"
set xlabel "t"
set ylabel "Cf"
plot file using 1:5 with lines linewidth 2 title "Cf"

unset multiplot
pause -1 "Press Enter to exit"