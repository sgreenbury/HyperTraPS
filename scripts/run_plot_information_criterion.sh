#!/bin/bash

prefix=""
pyscript=../src/python/plot_information_criterion.py
${pyscript}\
    -f "${prefix}forwards_information-criterion-0.csv"\
    -fontsize 6\
    -smooth_so 0\
    -outfile "${prefix}forwards-information-critertion"

# ../python/plot_information_criterion.py\
#     -f "${prefix}forwards_so_information-criterion-0.csv"\
#     -f2 "${prefix}forwards_fo_information-criterion-0.csv"\
#     -fontsize 6\
#     -smooth_so 0\
#     -outfile "${prefix}forwards-information-criterion-2"

