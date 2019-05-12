#!/bin/bash

prefix=""
pyscript=../src/python/plot_ordering_histogram.py
${pyscript}\
    -f "${prefix}forwards_list-pord-zero-one.csv"\
    -transition_data "${prefix}transitions.txt"\
    -labels "labels.csv"\
    -fontsize 6\
    -xevery 1\
    -outfile "forwards-ws1"

${pyscript}\
    -f "${prefix}forwards_list-pord-zero-one.csv"\
    -f2 "${prefix}forwards_list-pord-match-data.csv"\
    -transition_data "${prefix}transitions.txt"\
    -labels "labels.csv"\
    -fontsize 6\
    -xevery 1\
    -outfile "forwards-ws1-ws2"
