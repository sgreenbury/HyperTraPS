#!/bin/bash

../bin/RUN_PW -h

time ../bin/RUN_PW\
     -w zero-one\
     -b 50000\
     -R 100

time ../bin/RUN_PW\
     -w match-data\
     -b 50000\
     -R 100
