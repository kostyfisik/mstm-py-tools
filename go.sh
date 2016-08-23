#!/bin/bash
cd src
./mstm-input.py
cp mstm.inp ../bin
cd ../bin
rm *.dat
time mpirun -np 2 ./mstm mstm.inp
cp ../src/plot-field.py ./
./plot-field.py
#mpirun -np 2 ./mstm single-sphere.inp
#cat nf*
