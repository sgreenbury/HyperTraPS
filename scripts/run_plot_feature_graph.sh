#!/bin/bash

pyscript=../src/python/plot_feature_graph.py
prefix=""
# ${pyscript}\
#     -f "forwards.txt"\
#     -prob_type "joint"\
#     -layout_type "circular"\
#     -data_type "match-data"\
#     -width 3.4


${pyscript}\
    -f "forwards.txt"\
    -prob_type "joint"\
    -layout_type "circular"\
    -data_type "zero-one"\
    -width 4\
    -any_time 0\
    -node_size 100\
    -connection_style "arc3,rad=-0.3"\

# ${pyscript}\
#     -f "forwards.txt"\
#     -prob_type "joint"\
#     -layout_type "circular"\
#     -data_type "zero-one"\
#     -any_time 1\
#     -node_size 100\
#     -connection_style "arc3,rad=-0.3"\
#     -width 3.4

# ${pyscript}\
#     -f "forwards.txt"\
#     -prob_type "conditional"\
#     -layout_type "circular"\
#     -data_type "zero-one"\
#     -width 3.4
