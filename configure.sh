#!/bin/bash

# Set-up object and binary folders
rm -rf bin/
rm -rf obj/

mkdir bin/
mkdir obj/

# Compile C++ src
make all
