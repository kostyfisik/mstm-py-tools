#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import mstm_input as mi
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import os
from subprocess import call
def GetIndex(WLs, fname):
    data = np.loadtxt(fname)
    WL = data[:,0]*1000.0
    indexRe = data[:,1]
    indexIm = np.zeros(data.shape[0])+1e-4
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
    mstm_input.spheres.Reset()
    if R1 != 0:
        mstm_input.spheres.AddSphere(R1, [0, 0, 0], n1)
    if R2 != 0:
        mstm_input.spheres.AddSphere(R2, [R1+R2+Sep, 0, 0], n2)
###############################################################################
def ReadData(in_data,mstm_input, span):
    out_data = in_data
    isData = False
    WL = mstm_input.WL
    with open(mstm_input.sign+'.dat', 'r') as data_file:
        for data_line in data_file:
            # print(data_line)
            if "Qsca" in data_line:
                if len(out_data) == 0:
	            out_data += "# x "+data_line
	        isData = True
	        continue
	    if isData == True:
		if "0" in data_line:
		    out_data += str(span)+" "+data_line
		else:
		    isData = False
    return out_data
###############################################################################
def LoadData(fname):
    data = np.loadtxt(fname)
    number_of_spheres = int(max(data[:,1]))
    WLs = len(data[:,0])/number_of_spheres
    items = len(data[0,:])
    data_spaced = np.empty([number_of_spheres, WLs, items])
    # print(data_spaced)
    for j in range(len(data[:,0])):        
        row = data[j,:]
        data_spaced[ int(row[1])-1, int(j/number_of_spheres), :] = row
    return data, data_spaced
###############################################################################
def SaveSpectra(fname, WL, span, total_points):
    WLs = np.linspace(WL, WL, 1)
    index_BaTiO3 = GetIndex(WLs, "BaTiO3-Wemple-o.txt")[0][1]
    index_Au = GetIndex(WLs, "Au-Jhonson.txt")[0][1]
    #index_Au = GetIndex(WLs, "Au-Rakic.txt")    
    out_data=""
    mstm_input = mi.InputFile()
    R1,R2 = 0,0
    rg = np.linspace(-span, span, total_points)

    #R1 = 240/2 +300
    R1 = 407.5
    #R2 = 200/2
    Sep = 0
    for i in range(len(rg)):
        n1 = index_BaTiO3
        n2 = index_Au
        SetMstmModel(mstm_input, WL, R1+rg[i], R2, n1, n2, Sep)
        mstm_input.WriteFile()
        sign = mstm_input.sign
        print(sign)
        
        with  open(os.devnull, 'w') as FNULL:
            call(["mpirun", "-np", "2", "./mstm", "mstm.inp"],
                     stdout=FNULL)
            out_data = ReadData(out_data, mstm_input, R1+rg[i])
    with open(fname, 'w') as f:
			f.write(out_data)
    return sign
            


fname="spectra.dat"
WL = 600
span = 0.5 #nm
total_points = 201

sign = SaveSpectra(fname, WL, span, total_points)
sign = sign[:-6]+str(WL)+"nm"
data_all, data = LoadData(fname)
    # mstm_input.cut_plane = 'xy'
    # mstm_input.plot_scale = (R1+3*Sep+2*R2)/R1
    # mstm_input.points = 150
    # #mstm_input.isPlotField = False
    # mstm_input.isPlotField = True

plots=len(data)
############################# Plotting ######################
fig, axs = plt.subplots(plots,figsize=(4,2*plots), sharex=True)#, sharey=True)
# for ax in axs:
#     ax.locator_params(axis='y',nbins=4)
    # for label in ['left', 'right', 'top', 'bottom']:
    #     ax.spines[label].set_position(('outward',-1.3))
    #ax.tick_params(axis='x', pad=30)
for i in range(plots):
    plotwidth=.5
    ax = axs
    if plots > 1: ax = axs[i]
    cax = ax.plot(data[i,:,0], data[i,:,10], linewidth=plotwidth,
                         solid_joinstyle='round', solid_capstyle='round', color='black'
                         , label=r"$Q_{sca}$"
    )
    cax = ax.plot(data[i,:,0], data[i,:,11]/max(data[i,:,11])*max(data[i,:,9])/data[i,:,0]*max(data[i,:,0]), linewidth=plotwidth,
                         solid_joinstyle='round', solid_capstyle='round', color='red'
                         , label=r"$norm(Q_{abs})*max(Q_{ext})*max(R)/R$"
    )
    cax = ax.plot(data[i,:,0], data[i,:,9], linewidth=plotwidth/2,
                         solid_joinstyle='round', solid_capstyle='round', color='blue'
                         , label=r"$Q_{ext}$"
    )
    lg=ax.legend(loc='upper right',prop={'size':6})
    lg.draw_frame(False)

fname=sign
print(fname)
plt.savefig(fname+".pdf",pad_inches=0.02, bbox_inches='tight')

    # cax = axs[AgSi].plot(data[:,0], data[:,2], linewidth=plotwidth/1.5,
    #                      solid_joinstyle='round', solid_capstyle='round', color='red'
    #                      , label=r"$\tilde{b}_1$"
    # )
    # cax = axs[AgSi].plot(data[:,0], data[:,3], linewidth=plotwidth,
    #                      solid_joinstyle='round', solid_capstyle='round', color='green'
    #                      , label=r"$\tilde{a}_2$"
    # )
    # cax = axs[AgSi].plot(data[:,0], data[:,4], linewidth=plotwidth/1.5,
    #                      solid_joinstyle='round', solid_capstyle='round', color='blue'
    #                      , label=r"$\tilde{b}_2$"
    # )
    # axs[AgSi].axhline(y=0.25, ls='--', dashes=[2,2], color='gray')
    # lg=axs[AgSi].legend(loc='center right',prop={'size':11})
    # #lg=axs[SiAgSi].legend(loc='upper right',prop={'size':8})
    # #lg.get_frame().set_linewidth(0.0)
    # axs[AgSi].annotate('0.25', xy=(530, 0.25), fontsize=9, color='gray',
    #                 horizontalalignment='left', verticalalignment='bottom')

    # lg.draw_frame(False)

    # plotwidth=2.0
    # cax = axs[SiAgSi].plot(data_spaced2[:,0], data_spaced2[:,1], linewidth=plotwidth,
    #                      solid_joinstyle='round', solid_capstyle='round', color='black'
    #                      , label=r"$\tilde{a}_1$"
    # )

    # cax = axs[SiAgSi].plot(data_spaced2[:,0], data_spaced2[:,2], linewidth=plotwidth/1.5,
    #                      solid_joinstyle='round', solid_capstyle='round', color='red'
    #                      , label=r"$\tilde{b}_1$"
    # )
    # cax = axs[SiAgSi].plot(data_spaced2[:,0], data_spaced2[:,3], linewidth=plotwidth,
    #                      solid_joinstyle='round', solid_capstyle='round', color='green'
    #                      , label=r"$\tilde{a}_2$"
    # )
    # cax = axs[SiAgSi].plot(data_spaced2[:,0], data_spaced2[:,4], linewidth=plotwidth/1.5,
    #                      solid_joinstyle='round', solid_capstyle='round', color='blue'
    #                      , label=r"$\tilde{b}_2$"
    # )
    # axs[SiAgSi].axhline(y=0.25, ls='--', dashes=[2,2], color='gray')
    # lg=axs[SiAgSi].legend(loc='center right',prop={'size':11})
    # #lg=axs[SiSiAgSi].legend(loc='upper right',prop={'size':8})
    # #lg.get_frame().set_linewidth(0.0)
    # axs[SiAgSi].annotate('0.25', xy=(530, 0.25), fontsize=9, color='gray',
    #                 horizontalalignment='left', verticalalignment='bottom')

    # lg.draw_frame(False)

    # if isAll:
    #     plotwidth=2.0
    #     cax = axs[SiAgSi2].plot(data_spaced3[:,0], data_spaced3[:,1], linewidth=plotwidth,
    #                          solid_joinstyle='round', solid_capstyle='round', color='black'
    #                          , label=r"$\tilde{a}_1$"
    #     )

    #     cax = axs[SiAgSi2].plot(data_spaced3[:,0], data_spaced3[:,2], linewidth=plotwidth/1.5,
    #                          solid_joinstyle='round', solid_capstyle='round', color='red'
    #                          , label=r"$\tilde{b}_1$"
    #     )
    #     cax = axs[SiAgSi2].plot(data_spaced3[:,0], data_spaced3[:,3], linewidth=plotwidth,
    #                          solid_joinstyle='round', solid_capstyle='round', color='green'
    #                          , label=r"$\tilde{a}_2$"
    #     )
    #     cax = axs[SiAgSi2].plot(data_spaced3[:,0], data_spaced3[:,4], linewidth=plotwidth/1.5,
    #                          solid_joinstyle='round', solid_capstyle='round', color='blue'
    #                          , label=r"$\tilde{b}_2$"
    #     )
    #     axs[SiAgSi2].axhline(y=0.25, ls='--', dashes=[2,2], color='gray')
    #     lg=axs[SiAgSi2].legend(loc='center right',prop={'size':11})
    #     #lg=axs[SiSiAgSi2].legend(loc='upper right',prop={'size':8})
    #     #lg.get_frame().set_linewidth(0.0)
    #     axs[SiAgSi2].annotate('0.25', xy=(530, 0.25), fontsize=9, color='gray',
    #                     horizontalalignment='left', verticalalignment='bottom')

    #     lg.draw_frame(False)

    #     plotwidth=2.0
    #     # l,n = 0, 0
    #     # cax = axs[SiAgSi3].plot(data4[:,0], data4[:,l*6+n*2+1], linewidth=plotwidth/1.5,
    #     #                solid_joinstyle='round', solid_capstyle='round', color='red'
    #     #                , label="внутренний"
    #     #                )

    #     # l,n = 1, 0
    #     # cax = axs[SiAgSi3].plot(data4[:,0], data4[:,l*6+n*2+1], linewidth=plotwidth,
    #     #                solid_joinstyle='round', solid_capstyle='round', color='green'
    #     #                , label="средний"
    #     #                )
    #     # l,n = 2, 0
    #     # cax = axs[SiAgSi3].plot(data4[:,0], data4[:,l*6+n*2+1], linewidth=plotwidth,
    #     #                solid_joinstyle='round', solid_capstyle='round', color='black'
    #     #                , label='внешний'
    #     #                )


    #     # lg=axs[SiAgSi3].legend(loc='upper left',prop={'size':11})
    #     # lg.draw_frame(False)


    # y_up_lim = 0.29
    # axs[AgSi].set_ylabel(r'$\tilde{a}_n ,\ \tilde{b}_n$', labelpad=4.1)
    # axs[AgSi].set_ylim(0, y_up_lim)

    # axs[SiAgSi].set_ylabel(r'$\tilde{a}_n ,\ \tilde{b}_n$', labelpad=4.1)
    # axs[SiAgSi].set_ylim(0, y_up_lim)

    # axs[SiAgSi2].set_ylabel(r'$\tilde{a}_n ,\ \tilde{b}_n$', labelpad=4.1)
    # axs[SiAgSi2].set_ylim(0, y_up_lim)

    # # axs[SiAgSi3].set_ylabel(r'$\left|a_{ln}\right|^2+\left|d_{ln}\right|^2$', labelpad=+3)
    # axs[SiAgSi2].set_xlabel('Длина волны, нм', labelpad=2)
    # # axs[SiAgSi3].set_ylim(0.8, 8000)

    # plt.xlim(from_WL,  to_WL)

    # # axs[AgSi].annotate(r'$\Delta=%i$'%extra_width, xy=(0.09, 0.985), xycoords='axes fraction', fontsize=10,
    # #                 horizontalalignment='left', verticalalignment='top')

    # axs[AgSi].annotate('(г)', xy=(0.985, 0.98), xycoords='axes fraction', fontsize=10,
    #                 horizontalalignment='right', verticalalignment='top')
    # axs[SiAgSi].annotate('(д)', xy=(0.985, 0.98), xycoords='axes fraction', fontsize=10,
    #                 horizontalalignment='right', verticalalignment='top')
    # if isAll:
    #     axs[SiAgSi2].annotate('(е)', xy=(0.985, 0.98), xycoords='axes fraction', fontsize=10,
    #                 horizontalalignment='right', verticalalignment='top')
    #     # axs[SiAgSi3].annotate('(d)', xy=(0.99, 0.985), xycoords='axes fraction', fontsize=10,
    #     #             horizontalalignment='right', verticalalignment='top')
    #     # axs[SiAgSi3].set_yscale('log')
    #     # axs[SiAgSi3].set_yticks([1, 10, 100, 1000])
        
    #     # axs[SiAgSi3].tick_params(axis='y', which='major', pad=-0.65)

    #    # fmt = tkr.ScalarFormatter()
    #     # fmt.set_scientific(True)
    #     # axs[SiAgSi3].get_yaxis().set_major_formatter(fmt)

    # fig.subplots_adjust(hspace=.05)
    # plt.minorticks_off()
     # fname="2015-04-01-SiAgSi-ab-spectra4"
     #  plt.savefig(fname+".pdf",pad_inches=0.02, bbox_inches='tight')
    # #plt.draw()

    # #plt.show()

    # plt.clf()
    # plt.close()
