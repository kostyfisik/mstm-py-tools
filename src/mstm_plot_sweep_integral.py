#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import mstm_input as mi
import math
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.ticker as tkr 
###############################################################################
def LoadData(fname):
    data = np.loadtxt(fname)
    spheres_col_num = 2
    number_of_spheres = int(max(data[:,spheres_col_num]))
    WLs = len(data[:,0])/number_of_spheres
    items = len(data[0,:])
    data_spaced = np.empty([number_of_spheres, WLs, items])
    for j in range(len(data[:,0])):        
        row = data[j,:]
        data_spaced[ int(row[spheres_col_num])-1, int(j/number_of_spheres), :] = row
    return data, data_spaced
###############################################################################
def SaveSpectra(fname, WL,  from_span, to_span, total_points, b_ang, pol1, pol2):
    WLs = np.linspace(WL, WL/2.0, 2)
    out_data=""
    mstm_input = mi.InputFile()
    R1,R2 = 0,0
    rg = np.linspace( from_span, to_span, total_points)
    #mstm_input.isPlotField = True
    mstm_input.isPlotField = False
    #mstm_input.AddSource(WL, 0, 0, 0)
    #mstm_input.AddSource(WL, 0, 180, 0)
    
    R1 = 240/2 
    # R1 = 407.5965
    R1 = 580.5342
    #R1 = 350
    #R2 = 200/2
    #R2 = 100
    Sep = 10
    axis = 'z-'
    # n1 = "BaTiO3-Wemple-o"
    n1 = "BaTiO3-fixed"
    n2 = "Au-Jhonson"
    out_int = ""
    mstm_input.ResetSources()
    mstm_input.AddSource(WL, 0, 0, pol1)
    mstm_input.AddSource(WL, 0, b_ang, pol2)
    print("Span "+str(rg))
    for i in range(len(rg)):
        span_value = rg[i]
        mstm_input.ResetMstmModel(WL, span_value, R2, n1, n2, Sep,axis) # do span over R1
        sign = mstm_input.sign
        print(sign)
        out_data, integral = mstm_input.IntegrateOverlapNF(out_data, span_value)
        # if "N1" not in sign:
        #     sign += "--" + str(axis)
        out_int += str(span_value)+" "+str(integral)+"\n"
        print (span_value, integral)
    print(out_int)
    with open(fname, 'w') as f:
            f.write(out_data)
    with open(fname[:-4]+"-int.dat", 'w') as f:
            f.write(out_int)
    return sign

def run_and_plot(b_ang, pol1, pol2):
    fname="spectra.dat"
    WL = 600
    from_span = 50
    # span = 0.5 #nm
    # from_span =  580.5342-span
    to_span =  200#+span
    total_points = 501
    total_points = 16
    #total_points = 61
    sign = SaveSpectra(fname, WL, from_span, to_span, total_points, b_ang, pol1, pol2)
    sign = sign+"-"+str(WL)+"nm"
    data_int = np.loadtxt(fname[:-4]+"-int.dat")
    data_all, data = LoadData(fname)
    plots=len(data)
    ############################# Plotting ######################
    fig, axs = plt.subplots(2*plots,figsize=(4,2*2*plots), sharex=True)#, sharey=True)
    for i in range(plots):
        plotwidth=.5
        ax = axs[i]
        Q_ext, Q_sca, Q_abs = 10, 11, 12
        cax = ax.plot(data[i,0::2,0], data[i,0::2,Q_sca], linewidth=plotwidth,
                      solid_joinstyle='round', solid_capstyle='round', color='black'
                      , label=r"$Q_{sca}$"
        )
        cax = ax.plot(data[i,1::2,0], data[i,1::2,Q_sca], linewidth=plotwidth,
                      solid_joinstyle='round', solid_capstyle='round', color='black'
                      , label=r"$Q_{sca}$"
        )
        # cax = ax.plot(data[i,:,0], data[i,:,Q_ext], linewidth=plotwidth/2,
        #               solid_joinstyle='round', solid_capstyle='round', color='blue'
        #               , label=r"$Q_{ext}$"
        # )
        # cax = ax.plot(data[i,:,0], data[i,:,Q_abs], linewidth=plotwidth,
        #               solid_joinstyle='round', solid_capstyle='round', color='red'
        #               , label=r"$Q_{abs}$")
        cax = ax.plot(data_int[:,0], abs(data_int[:,1])/max(abs(data_int[:,1]))*max(data[i,:,Q_ext]),
                      linewidth=plotwidth*2,
                      solid_joinstyle='round', solid_capstyle='round', color='green'
                      , label=r"$norm(\frac{1}{V}\int{\,E_{0}(\omega)\,E_{1}(\omega)\,dV})$"
        )
        lg=ax.legend(loc='upper left',prop={'size':6})
        lg.draw_frame(False)
        ax = axs[i+plots]
        cax = ax.plot(data_int[:,0], data_int[:,1], linewidth=plotwidth*2,
                      solid_joinstyle='round', solid_capstyle='round', color='green'
                      , label=r"$\frac{1}{V}\int{\,E_{0}(\omega)\,E_{1}(\omega)\,dV}$"
        )
        lg=ax.legend(loc='upper left',prop={'size':6})
        lg.draw_frame(False)

    fname=sign
    plt.savefig(fname+".pdf",pad_inches=0.02, bbox_inches='tight')




b_ang_set = [0]
# b_ang_set = [0,45,90,135,180]
# b_ang_set = [0,30,45,60,90,120,135,150,180]
pol_ang_set = [0]
for b_ang in b_ang_set:
    for pol1 in pol_ang_set:
        for pol2 in pol_ang_set:
            run_and_plot(b_ang,pol1,pol2)
