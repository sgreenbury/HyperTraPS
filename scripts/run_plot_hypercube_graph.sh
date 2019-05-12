#!/bin/bash

pyscript=../src/python/plot_hypercube_graph.py
prefix=""
${pyscript}\
    -f "${prefix}forwards_list-pord-zero-one.csv"\
    -extra_labels_routes "20,2"\
    -label_type "greedy_data"\
    -labels_fontsize 4\
    -transition_data "${prefix}transitions.txt"\
    -labels "${prefix}labels.csv"\
    -aspect 0.9\
    -width 3.4\
    -gamma 0\
    -xlabel ""\
    -extra "no"\
    -xevery 1\
    -outfile_graph "${prefix}forwards-hypercube-graph-zero-one-no01"

${pyscript}\
    -f "${prefix}forwards_list-pord-match-data.csv"\
    -extra_labels_routes "20,2"\
    -label_type "greedy_data"\
    -labels_fontsize 4\
    -transition_data "${prefix}transitions.txt"\
    -labels "${prefix}labels.csv"\
    -aspect 0.9\
    -width 3.4\
    -gamma 0\
    -xlabel ""\
    -extra "no"\
    -xevery 1\
    -outfile_graph "${prefix}forwards-hypercube-graph-match-data-no01"

${pyscript}\
    -f "${prefix}forwards_list-pord-zero-one.csv"\
    -write_transitions "yes"\
    -extra_labels_routes "7,2"\
    -layout_type "spring"\
    -gamma 0\
    -label_type "greedy_data"\
    -labels_fontsize 6\
    -transition_data "${prefix}transitions.txt"\
    -labels "${prefix}labels.csv"\
    -aspect 0.9\
    -width 3.4\
    -gamma 0\
    -xlabel ""\
    -extra "no"\
    -xevery 1\
    -outfile_graph "${prefix}forwards-hypercube-graph-zero-one-g0"

${pyscript}\
    -f "${prefix}forwards_list-pord-match-data.csv"\
    -write_transitions "yes"\
    -layout_type "spring"\
    -gamma 0\
    -extra_labels_routes "7,2"\
    -label_type "greedy_data"\
    -labels_fontsize 6\
    -transition_data "${prefix}transitions.txt"\
    -labels "${prefix}labels.csv"\
    -aspect 0.9\
    -width 3.4\
    -gamma 0\
    -xlabel ""\
    -extra "no"\
    -xevery 1\
    -outfile_graph "${prefix}forwards-hypercube-graph-match-data-g0"
