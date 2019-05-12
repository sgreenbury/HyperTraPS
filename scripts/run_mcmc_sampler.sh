#!/bin/bash

../bin/RUN_MCMC_SAMPLER -h

time ../bin/RUN_MCMC_SAMPLER\
     -f transitions.txt\
     -M second-order\
     -N 1000\
     -r 20\
     -n zero\
     -p 20\
     -k mcmc-apm\
     -s 0.05\
     -b 50000\
     -i 100\
     -S 0\
     -t 3\
     -q 0
