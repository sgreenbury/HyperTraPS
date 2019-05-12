#!/bin/bash

../bin/RUN_PW -h

time ../bin/RUN_PW\
     -L forwards_regularised.txt\
     -w zero-one\
     -b 0\
     -e 100\
     -i 1\
     -g 1\
     -R 100

time ../bin/RUN_PW\
     -L forwards_regularised.txt\
     -w match-data\
     -b 0\
     -e 100\
     -i 1\
     -g 1\
     -R 100
