#!/bin/bash
cd src
./mstm-input.py
cp mstm.inp ../bin
cd ../bin
rm *.dat
cp ../src/*.py ./
cp ../data/Au-Jhonson.txt ./
cp ../data/BaTiO3-Wemple-o.txt ./

./mstm_plot_spectra.py
# time mpirun -np 2 ./mstm mstm.inp
# time ./plot-field.py
# time ./check-field.py
#mpirun -np 2 ./mstm single-sphere.inp
#cat nf*
