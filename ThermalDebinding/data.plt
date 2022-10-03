#!/usr/bin/env -S gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,14'
    set output ARG1.".eps"
} else {
    set terminal qt size 1000, 400 font 'Helvetica,14'
}

file = "log"
# Generate 3 columns
data = sprintf("<(awk '/y =/{print $3, $6, $9}' %s)", file)

set multiplot layout 1, 2

set xlabel "Temperature"
#set format y "%g"

set ylabel "Mass fraction"
plot data u 2:1 w l lw 2 title "polymer"

set ylabel "Pressure (Pa)"
set log y
plot data u 2:3 w l lw 2 title "maximum monomer pressure"

