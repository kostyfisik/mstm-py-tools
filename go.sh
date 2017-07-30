#!/bin/bash
cd bin
rm *.dat *.pdf
cp ../mstm-v3/mstm ./
cp ../src/*.py ./
cp ../data/Au*.txt ./
cp ../data/BaTiO3*.txt ./
cp ../data/Si*.txt ./

# ./mstm_plot_spectra.py
#./mstm_plot_sweep.py
#./mstm_plot_sweep_integral.py
./mstm_plot_spectra_Si_dimer.py

#ls
# time mpirun -np 2 ./mstm mstm.inp

#time ./plot_field.py

# time ./check-field.py
#mpirun -np 2 ./mstm single-sphere.inp
#cat nf*
