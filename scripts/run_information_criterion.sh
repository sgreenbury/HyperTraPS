#!/bin/bash

../bin/RUN_INFORMATION_CRITERION -h

time ../bin/RUN_INFORMATION_CRITERION\
     -f transitions.txt\
     -t 3\
     -r 20\
     -L forwards.txt\
     -T 100\
     -G 0\
     -I AIC
