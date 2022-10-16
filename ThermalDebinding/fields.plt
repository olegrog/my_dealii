#!/usr/bin/env -S gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,14'
    set output ARG1.".eps"
} else {
    set terminal qt size 1000, 700 font 'Helvetica,14'
}

file = "log"
# Generate 3 columns
data1 = sprintf("<(awk '/T =/{print $3, $6, $9}' %s)", file)
data2 = sprintf("<(awk '/T =/{T=$3} /max/{print T, $3, $6, $9, $12, $15, $18}' %s)", file)

set multiplot layout 2, 2

set xlabel "Temperature (°C)"

KtoC(T) = T - 273.15

set ylabel "Mass fraction of a polymer specie"
set yrange [0:1]
plot data1 u (KtoC($1)):2 w l lw 2 title "Polymer1", \
  '' u (KtoC($1)):3 w l lw 2 title "Polymer2"

set ylabel "Mass diffusivity (m^2/s)"
set log y
set yrange [1e-15:1e-5]
set style fill transparent solid 0.5
plot data2 u (KtoC($1)):6:7 w filledcurves title "Due to convection (min/max)", \
    data2 u (KtoC($1)):4:5 w filledcurves title "Due to diffusion (min/max)", \
    data2 u (KtoC($1)):4 w l lw 2 lt -1 notitle, \
    data2 u (KtoC($1)):5 w l lw 2 lt -1 notitle, \
    data2 u (KtoC($1)):6 w l lw 2 lt -1 notitle, \
    data2 u (KtoC($1)):7 w l lw 2 lt -1 notitle


set ylabel "Maximum density of monomer (kg/m^3)"
unset yrange
plot data2 u (KtoC($1)):3 w l notitle lw 2

set ylabel "Pressure (Pa)"
plot data2 u (KtoC($1)):2 w l lw 2 title "Maximum monomer pressure"

