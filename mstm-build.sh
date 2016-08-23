#!/bin/bash
cd mstm-v3
pwd
mpif90 -g -o mstm mpidefs-parallel-v3.0.f90 mstm-intrinsics-v3.0.f90 mstm-modules-v3.0.f90 mstm-main-v3.0.f90
