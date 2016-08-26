#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import mstm_input as mi
import math
import numpy as np
def GetIndex(WLs, fname):
    data = np.loadtxt(fname)
    WL = data[:,0]*1000.0
    indexRe = data[:,1]
    indexIm = np.zeros(data.shape[0])
    if data.shape[1] == 3:
	indexIm = data[:,2]
    from scipy.interpolate import interp1d
    fRe = interp1d(WL, indexRe)
    fIm = interp1d(WL, indexIm)
    # fRe = interp1d(WL, indexRe, kind=2)
    # fIm = interp1d(WL, indexIm, kind=2)

    data = np.vstack((WLs, fRe(WLs)+fIm(WLs)*1j))
    # data = np.concatenate(WLs, np.array(fRe(WLs)))
    # data = np.concatenate(WLs, )
    return np.transpose(data)
###############################################################################
def SetMstmModel(mstm_input, WL, R1, R2, n1, n2, Sep):
    mstm_input.WL = WL
    mstm_input.spheres.AddSphere(R1, [0, 0, 0], n1)
    if R2 != 0:
        mstm_input.spheres.AddSphere(R2, [R1+R2+Sep, 0, 0], n2)
	
###############################################################################
def SaveSpectra(fname, from_WL, to_WL, total_points):
    WLs = np.linspace(from_WL, to_WL, total_points)
    index_BaTiO3 = GetIndex(WLs, "BaTiO3-Wemple-o.txt")
    index_Au = GetIndex(WLs, "Au-Jhonson.txt")    
    print(index_Au)
    WL=400
    R1 = 200
    n1 = 2+0.1j
    n2 = 3+2j
    R2 = 100
    Sep = 10

    mstm_input = mi.InputFile()
    SetMstmModel(mstm_input, WL, R1, R2, n1, n2, Sep):
    mstm_input.WriteFile()

    print(mstm_input.PrintInput())


fname="spectra.dat"
from_WL = 400
to_WL = 800
total_points = 3

SaveSpectra(fname, from_WL, to_WL, total_points)

    # mstm_input.cut_plane = 'xy'
    # mstm_input.plot_scale = (R1+3*Sep+2*R2)/R1
    # mstm_input.points = 150
    # #mstm_input.isPlotField = False
    # mstm_input.isPlotField = True
