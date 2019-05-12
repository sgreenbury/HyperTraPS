#!/bin/bash

prefix=""
pyscript=../src/python/plot_raw_data.py
$pyscript\
    -f raw-data.txt\
    -aspect 0.8\
    -fontsize 6\
    -labels labels.csv\
    -row_every 1\
    -outfile_type "pdf"\
    -xlabel "Samples"
