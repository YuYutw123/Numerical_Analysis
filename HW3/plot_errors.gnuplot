# Set output format and terminal
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output 'pC_3.png'

# Title and labels
set title "Error Norms" font ",14"
set xlabel "Polynomial Degree" font ",12"
set ylabel "Norm Value" font ",12"

# Grid and style
set grid
set key outside top right

# Line styles
set style line 1 lc rgb '#0072bd' lw 2 pt 7 ps 1.5 # Blue line for 2-norm
set style line 2 lc rgb '#d95319' lw 2 pt 5 ps 1.5 # Red line for infinity-norm

# Load data and plot
plot 'data.txt' using 1:2 with linespoints ls 1 title "2-Norm", \
     'data.txt' using 1:3 with linespoints ls 2 title "Infinity-Norm"
