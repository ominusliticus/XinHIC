#!/bin/bash

g++-11 -Wall -O3 -std=c++17 -o main.x MeatAndBones.cpp main.cpp /usr/local/lib/libfmt.a
./main.x

python3 make_plots.py
# evince trial_file.pdf 

cp D-self-energies.pdf ../Docs/Figures/
cp D-self-energies.pdf ../Docs/Figures/Figures
