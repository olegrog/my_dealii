#!/bin/bash -e

logfile=log
[[ $# -eq 1 ]] && logfile=$1

if [ -z $DISPLAY ]; then
    cmd='set term dumb;'
fi

gnuplot -p -e "$cmd set format y '%.1e'; set log y; plot '<&3' w lp" \
    3< <(grep delta_t $logfile | awk '{print $4}')
