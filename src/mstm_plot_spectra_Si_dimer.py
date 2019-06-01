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
    WLs = int(len(data[:,0])/number_of_spheres)+1
    items = len(data[0,:])
    print(number_of_spheres,
         WLs,
         items)
    data_spaced = np.empty(
        (number_of_spheres,
         WLs,
         items))
    for j in range(len(data[:,0])):        
        row = data[j,:]
        data_spaced[ int(row[spheres_col_num])-1, int(j/number_of_spheres), :] = row
    return data, data_spaced
###############################################################################
def SaveSpectra(fname, from_WL, to_WL, WL_points):
    WLs = np.linspace(from_WL, to_WL, WL_points)
    out_data=""
    mstm_input = mi.InputFile()
    R1,R2 = 0,0
    #mstm_input.isPlotField = True
    mstm_input.isPlotField = False
    #mstm_input.AddSource(WL, 0, 0, 0)
    #mstm_input.AddSource(WL, 0, 180, 0)
    
    R1 = 50
    R2 = 50
    Sep = 100
    axis = 'x'
    n1 = "Si-Belyakov"
    n2 = "Si-Belyakov"
    out_int = ""
    for WL in WLs:
        mstm_input.ResetSources()
        mstm_input.AddSource(WL, 0, 0, 0)
        mstm_input.ResetMstmModel(WL, R1, R2, n1, n2, Sep,axis) # do span over R1
        sign = mstm_input.sign
        print(sign)
        out_data = mstm_input.RunForQsca(out_data)
        #out_data, integral = mstm_input.IntegrateOverlapNF(out_data, span_value)
        # out_int += str(span_value)+"\n"
        # print(out_int)
    print (out_data)
    with open(fname, 'w') as f:
            f.write(out_data)
    with open(fname[:-4]+"-int.dat", 'w') as f:
            f.write(out_int)
    return sign

def run_and_plot():
    fname="spectra.dat"
    from_WL = 350
    to_WL = 550
    WL_points = 21
    from_span = 100
    # span = 0.5 #nm
    # from_span =  580.5342-span
    to_span =  100#+span
    total_points = 501
    total_points = 1
    #total_points = 61
    sign = SaveSpectra(fname,  from_WL, to_WL, WL_points)
    sign = sign+"-"+str(from_WL)+"-"+str(to_WL)+"nm"
    data_int = np.loadtxt(fname)
    data_all, data = LoadData(fname)
    plots=len(data)
    # ############################# Plotting ######################
    # fig, axs = plt.subplots(2*plots,figsize=(4,2*2*plots), sharex=True)#, sharey=True)
    # for i in range(plots):
    #     plotwidth=.5
    #     ax = axs[i]
    #     Q_ext, Q_sca, Q_abs = 10, 11, 12
    #     cax = ax.plot(data[i,0::2,0], data[i,0::2,Q_sca], linewidth=plotwidth,
    #                   solid_joinstyle='round', solid_capstyle='round', color='black'
    #                   , label=r"$Q_{sca}$"
    #     )
    #     cax = ax.plot(data[i,1::2,0], data[i,1::2,Q_sca], linewidth=plotwidth,
    #                   solid_joinstyle='round', solid_capstyle='round', color='black'
    #                   , label=r"$Q_{sca}$"
    #     )
    #     # cax = ax.plot(data[i,:,0], data[i,:,Q_ext], linewidth=plotwidth/2,
    #     #               solid_joinstyle='round', solid_capstyle='round', color='blue'
    #     #               , label=r"$Q_{ext}$"
    #     # )
    #     # cax = ax.plot(data[i,:,0], data[i,:,Q_abs], linewidth=plotwidth,
    #     #               solid_joinstyle='round', solid_capstyle='round', color='red'
    #     #               , label=r"$Q_{abs}$")
    #     cax = ax.plot(data_int[:,0], abs(data_int[:,1])/max(abs(data_int[:,1]))*max(data[i,:,Q_ext]),
    #                   linewidth=plotwidth*2,
    #                   solid_joinstyle='round', solid_capstyle='round', color='green'
    #                   , label=r"$norm(\frac{1}{V}\int{\,E_{0}(\omega)\,E_{1}(\omega)\,dV})$"
    #     )
    #     lg=ax.legend(loc='upper left',prop={'size':6})
    #     lg.draw_frame(False)
    #     ax = axs[i+plots]
    #     cax = ax.plot(data_int[:,0], data_int[:,1], linewidth=plotwidth*2,
    #                   solid_joinstyle='round', solid_capstyle='round', color='green'
    #                   , label=r"$\frac{1}{V}\int{\,E_{0}(\omega)\,E_{1}(\omega)\,dV}$"
    #     )
    #     lg=ax.legend(loc='upper left',prop={'size':6})
    #     lg.draw_frame(False)

    # fname=sign
    # plt.savefig(fname+".pdf",pad_inches=0.02, bbox_inches='tight')

##############################################################################
# end of
# def run_and_plot(b_ang, pol1, pol2):
##############################################################################



# b_ang_set = [0]
# # b_ang_set = [0,45,90,135,180]
# # b_ang_set = [0,30,45,60,90,120,135,150,180]
# pol_ang_set = [0]
# for b_ang in b_ang_set:
#     for pol1 in pol_ang_set:
#         for pol2 in pol_ang_set:
#             run_and_plot(b_ang,pol1,pol2)
run_and_plot()
