#!/usr/bin/env -S gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,14'
    set output ARG1.".eps"
} else {
    set terminal qt size 1000, 400 font 'Helvetica,14'
}

file = "log"
# Generate 3 columns
data = sprintf("<(awk '/T =/{print $3, $6, $9, $12}' %s)", file)

set multiplot layout 1, 2

set xlabel "Temperature (°C)"

KtoC(T) = T - 273.15

#set format y "%g"

set ylabel "Mass fraction of a polymer specie"
set yrange [0:1]
plot data u (KtoC($1)):3 w l lw 2 title "Polymer1", '' u (KtoC($1)):4 w l lw 2 title "Polymer2"

set ylabel "Pressure (Pa)"
unset yrange
set log y
plot data u (KtoC($1)):2 w l lw 2 title "Maximum monomer pressure"

