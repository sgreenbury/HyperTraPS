#!/bin/bash

pyscript=../src/python/convert_data_to_transitions.py
if [ "$#" -eq 2 ]; then
    ${pyscript}\
	-data $1\
	-input_type $2
fi

if [ "$#" -eq 3 ]; then
    ${pyscript}\
	-data $1\
	-input_type $2\
	-phylogeny $3
fi
