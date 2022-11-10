#!/usr/bin/env -S gnuplot --persist -c

if (ARGC > 0) {
    set terminal postscript eps color font 'Helvetica,11'
    set output ARG1.".eps"
    label_shift = 0.05
} else {
    label_shift = 0
    set terminal qt size 1000, 700 font 'Helvetica,14'
}

file = "log"
# Generate 3 columns
data1 = sprintf("<(awk '/^T =/{print $3, $6, $9, $12, $15}' %s)", file)
data2 = sprintf("<(awk '/^T =/{T=$3} /^max/{print T, $3, $6, $9, $12, $15, $18}' %s)", file)
data3 = sprintf("<(awk '/^T =/{T=$3} /^organic/{print T, $3}' %s)", file)
data4 = "strength.txt"

set encoding utf8
set multiplot layout 2, 2
set xlabel "Temperature (°C)"

KtoC(T) = T - 273.15
fromMPa(p) = p*1e6
filter(x,min,max) = (x > min && x < max) ? x : 1/0

unset label; set label "a)" at graph -0.11 - label_shift, 1
set ylabel "Mass and volume fractions"
set yrange [0:1]
plot data3 u (KtoC($1)):2 w l title "Organic relative mass" lw 2, \
  data1 u (KtoC($1)):3 w l lw 2 title "Polymer1", \
  '' u (KtoC($1)):4 w l lw 2 title "Polymer2", \
  '' u (KtoC($1)):($3+$4) w l title "All polymers" lw 2, \
  '' u (KtoC($1)):2 w l lw 2 title "Porosity"

unset label; set label "b)" at graph -0.15 - label_shift, 1
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


unset label; set label "c)" at graph -0.12 - label_shift, 1
set ylabel "Mass concentration (kg/m^3)"
set yrange [1e-2:1e3]
plot data2 u (KtoC($1)):3 w l title "Maximum of monomer" lw 2, \
    data1 u (KtoC($1)):5 w l title "Polymer" lw 2

unset label; set label "d)" at graph -0.15 - label_shift, 1
set ylabel "Pressure (Pa)"
set yrange [1e4:1e9]
set style fill transparent solid 1
plot data2 u (filter(KtoC($1), 100, 200)):2 with filledcurves x1 notitle lc rgb "gray",\
    data2 u (KtoC($1)):2 w l title "Maximum of monomer pressure" lw 2, \
    data4 u 1:(fromMPa($2)) w lp title "Compressive strength", \
    data4 u 1:(fromMPa($3)) w lp title "Tensile strength", \
    data4 u 1:(fromMPa($4)) w lp title "Flexural strength", \

getTemperature(row) = system('awk ''/^T =/{print $3}'' '.file.' | sed '''.row.'q;d''')
stats data2 u 2 nooutput
print "Max pressure (MPa) = ", STATS_max/1e6, " when T (°C) = ", KtoC(getTemperature(STATS_index_max + 1))

