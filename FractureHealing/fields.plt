#!/usr/bin/env -S gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,11'
    set output ARG1.".eps"
    label_shift = 0.05
} else {
    label_shift = 0
    set terminal qt size 1000, 700 font 'Helvetica,14'
}

file1 = 'integral_values.dat'

set encoding utf8
set multiplot layout 3, 1
set xlabel "Time"
set key autotitle columnheader
set datafile separator tab
set datafile commentschars ""

set ylabel "Number of cells"
plot for [i=2:4] file1 u 1:i w l lw 2 title word(columnheader(i), 1)

set ylabel "Extracellular matrix mass"
plot for [i=5:6] file1 u 1:i w l lw 2 title word(columnheader(i), 1)

set ylabel "Growth factor mass"
plot for [i=7:8] file1 u 1:i w l lw 2 title word(columnheader(i), 1)

