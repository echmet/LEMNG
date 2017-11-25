#!/usr/bin/gnuplot

set encoding utf8
set decimalsign locale "en_US.UTF-8"
set terminal pngcairo size 3200,3200
set output "efg_plot.png"

set datafile separator ";"

plot [:][:] "efgplots.csv" using 1:2 with lines, \
            "efgplots.csv" using 1:3 with lines, \
            "efgplots.csv" using 1:4 with lines, \
            "efgplots.csv" using 1:5 with lines, \
            "efgplots.csv" using 1:6 with lines, \
            "efgplots.csv" using 1:8 with lines, \
            "efgplots.csv" using 1:9 with lines,

