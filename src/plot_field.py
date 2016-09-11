#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2009-2015 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Several functions to plot field and streamlines (power flow lines).

import numpy as np
import cmath


def GetCSTField(data_txt, WL):
    skips = 0
    with open(data_txt, 'r') as data_file:
        for data_line in data_file:
            if len(data_line.split()) > 4: break
            skips += 1
    # print("Skips = "+ str(skips))
    x=np.transpose(np.loadtxt(data_txt, skiprows=skips))
    E = np.sqrt(np.absolute(x[2]+1.0j*x[3])**2+np.absolute(x[4]+1.0j*x[5])**2+np.absolute(x[6]+1.0j*x[7])**2)
    H = np.sqrt(np.absolute(x[8]+1.0j*x[9])**2+np.absolute(x[10]+1.0j*x[11])**2+np.absolute(x[12]+1.0j*x[13])**2)
    coordX = np.unique(x[0])/WL
    coordZ = np.unique(x[1])/WL
    # print(x[2])
    return E, H, coordX, coordZ

def fieldplot(fig, ax, WL, filetoplot, comment='', WL_units=' ', crossplane='XZ',
              field_to_plot='Pabs', npts=101, factor=2.1, flow_total=11,
              is_flow_extend=True, pl=-1, outline_width=1, subplot_label=' '):
    #cst_data_txt="mie-fine.txt"
    #cst_data_txt="mie.txt"
    
    cst_data_txt=filetoplot #"N2__R_200_100__D_310_Sep_10__nf.dat"
    Er, Hr, coordZ, coordX = GetCSTField(cst_data_txt,WL)
    # print(coordZ)
    npts = len(coordX)
    try:
        from matplotlib import cm
        from matplotlib.colors import LogNorm

        if field_to_plot == 'Eabs':
            Eabs = Er
            Eabs_data = np.fliplr(np.resize(Eabs, (len(coordX), len(coordZ))).T)
            label = r'$|E|$'
        if field_to_plot == 'Habs':
            Eabs = Hr
            Eabs_data = np.fliplr(np.resize(Eabs, (len(coordX), len(coordZ))).T)
            label = r'$|H|$'
        # # WL units
        # WL=1
        # WL_units = r'$\lambda$'
        # scale_x = np.linspace(
        #     (min(coordX)), (max(coordX)), len(Eabs_data[0]))
        # scale_z = np.linspace(min(coordZ), max(coordZ), len(coordZ))
        # ax.locator_params(nbins=5)

        # # Define scale ticks
        # min_tick = np.amin(Eabs_data[~np.isnan(Eabs_data)])
        # max_tick = np.amax(Eabs_data[~np.isnan(Eabs_data)])
        # # min_tick = 0.0
        # # max_tick = 2.59
        # scale_ticks = np.linspace(min_tick, max_tick, 10)
        # #scale_ticks = np.power(10.0, np.linspace(np.log10(min_tick), np.log10(max_tick), 6))
        # #scale_ticks = [0.1,0.3,1,3,10, max_tick]
        # # Interpolation can be 'nearest', 'bilinear' or 'bicubic'
        # # ax.set_title(label)
        # # build a rectangle in axes coords
        # ax.annotate(subplot_label, xy=(0.0, 1.1), xycoords='axes fraction',  # fontsize=10,
        #             horizontalalignment='left', verticalalignment='top')
        # ax.xaxis.set_tick_params(width=outline_width/2.0)
        # ax.yaxis.set_tick_params(width=outline_width/2.0)
        # # ax.text(right, top, subplot_label,
        # #         horizontalalignment='right',
        # #         verticalalignment='bottom',
        # #         transform=ax.transAxes)
        cax = ax.imshow(Eabs_data, interpolation='none', cmap=cm.rainbow )

        # cax = ax.imshow(Eabs_data, interpolation='none', cmap=cm.rainbow,
        #                 origin='lower', vmin=min_tick, vmax=max_tick, extent=(min(scale_x), max(scale_x), min(scale_z), max(scale_z))
        #                 # ,norm = LogNorm()
        #                 )
        # ax.axis("image")

        # # Add colorbar
        # cbar = fig.colorbar(cax, ticks=[a for a in scale_ticks], ax=ax)
        # # vertically oriented colorbar
        # if 'angle' in field_to_plot:
        #     cbar.ax.set_yticklabels(['%3.0f' % (a) for a in scale_ticks])
        # else:
        #     cbar.ax.set_yticklabels(['%3.2f' % (a) for a in scale_ticks])
        # # pos = list(cbar.ax.get_position().bounds)
        # #fig.text(pos[0] - 0.02, 0.925, '|E|/|E$_0$|', fontsize = 14)
        # lp2 = -5.0
        # lp1 = -1.0
        # if crossplane == 'XZ':
        #     ax.set_xlabel('Z, ' + WL_units, labelpad=lp1)
        #     ax.set_ylabel('X, ' + WL_units, labelpad=lp2)
        # elif crossplane == 'YZ':
        #     ax.set_xlabel('Z, ' + WL_units, labelpad=lp1)
        #     ax.set_ylabel('Y, ' + WL_units, labelpad=lp2)
        # elif crossplane=='XYZ':
        #     ax.set_xlabel(r'$Z,\lambda$'+WL_units, fontsize = 25)
        #     ax.set_ylabel(r'$Y:X,\lambda$'+WL_units, fontsize = 25, labelpad=lp2)
        #     ax.axhline(y=0.0, ls='--', dashes=[7,5], color='gray', lw=outline_width)
        #     bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="white", lw=2)
        #     ax.annotate('E-k', xy=(0.95, 0.95), xycoords='axes fraction', fontsize=21,
        #                 horizontalalignment='right', verticalalignment='top',bbox = bbox_props)
        #     ax.annotate('H-k', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=21,
        #                 horizontalalignment='right', verticalalignment='bottom',bbox = bbox_props)

        # elif crossplane == 'XY':
        #     ax.set_xlabel('Y, ' + WL_units, labelpad=lp1)
        #     ax.set_ylabel('X, ' + WL_units, labelpad=lp2)
        # # # This part draws the nanoshell
        # from matplotlib import patches
        # from matplotlib.path import Path
        # for xx in [1.5/2.0, 1.5/2.0+0.8/3.75]:
        # #for xx in x:
        #     r = xx 
        #     s1 = patches.Arc((0, 0), 2.0 * r, 2.0 * r,  angle=0.0, zorder=1.8,
        #                      theta1=0.0, theta2=360.0, linewidth=outline_width, color='black')
        #     ax.add_patch(s1)
#         #
#         # for flow in range(0,flow_total):
#         #     flow_x, flow_z = GetFlow(scale_x, scale_z, Ec, Hc,
#         #                              min(scale_x)+flow*(scale_x[-1]-scale_x[0])/(flow_total-1),
#         #                              min(scale_z),
#         #                              #0.0,
#         #                              npts*16)
#         #     verts = np.vstack((flow_z, flow_x)).transpose().tolist()
#         #     #codes = [Path.CURVE4]*len(verts)
#         #     codes = [Path.LINETO]*len(verts)
#         #     codes[0] = Path.MOVETO
#         #     path = Path(verts, codes)
#         #     patch = patches.PathPatch(path, facecolor='none', lw=1, edgecolor='yellow')
#         #     ax.add_patch(patch)
#         if (not crossplane == 'XY') and flow_total > 0:

#             from matplotlib.path import Path
#             scanSP = np.linspace(-factor * x[-1], factor * x[-1], npts)
#             min_SP = -factor * x[-1]
#             step_SP = 2.0 * factor * x[-1] / (flow_total - 1)
#             x0, y0, z0 = 0, 0, 0
#             max_length = factor * x[-1] * 10
#             # max_length=factor*x[-1]*5
#             max_angle = np.pi / 160
#             if is_flow_extend:
#                 rg = range(0, flow_total * 5 + 1)
#             else:
#                 rg = range(0, flow_total)
#             for flow in rg:
#                 if is_flow_extend:
#                     f = min_SP*2 + flow*step_SP
#                 else:
#                     f = min_SP + flow*step_SP
#                 if crossplane=='XZ':
#                     x0 = f
#                 elif crossplane=='YZ':
#                     y0 = f
#                 elif crossplane=='XYZ':
#                     x0 = 0
#                     y0 = 0
#                     if f > 0:
#                         x0 = f
#                     else:
#                         y0 = f
#                 z0 = min_SP
#                     # x0 = x[-1]/20
#                 flow_xSP, flow_ySP, flow_zSP = GetFlow3D(
#                     x0, y0, z0, max_length, max_angle, x, m, pl)
#                 if crossplane == 'XZ':
#                     flow_z_plot = flow_zSP * WL / 2.0 / np.pi
#                     flow_f_plot = flow_xSP * WL / 2.0 / np.pi
#                 elif crossplane == 'YZ':
#                     flow_z_plot = flow_zSP * WL / 2.0 / np.pi
#                     flow_f_plot = flow_ySP * WL / 2.0 / np.pi
#                 elif crossplane=='XYZ':
#                     if f > 0:
#                         flow_z_plot = flow_zSP*WL/2.0/np.pi
#                         flow_f_plot = flow_xSP*WL/2.0/np.pi
#                     else:
#                         flow_z_plot = flow_zSP*WL/2.0/np.pi
#                         flow_f_plot = flow_ySP*WL/2.0/np.pi

#                 verts = np.vstack(
#                     (flow_z_plot, flow_f_plot)).transpose().tolist()
#                 codes = [Path.LINETO] * len(verts)
#                 codes[0] = Path.MOVETO
#                 path = Path(verts, codes)
#                 #patch = patches.PathPatch(path, facecolor='none', lw=0.2, edgecolor='white',zorder = 2.7)
#                 patch = patches.PathPatch(
#                     path, facecolor='none', lw=outline_width, edgecolor='white', zorder=1.9, alpha=0.7)
#                 # patch = patches.PathPatch(
#                 #     path, facecolor='none', lw=0.7, edgecolor='white', zorder=1.9, alpha=0.7)
#                 ax.add_patch(patch)
# #                ax.plot(flow_z_plot, flow_f_plot, 'x', ms=2, mew=0.1,
# #                        linewidth=0.5, color='k', fillstyle='none')

    finally:
        # terms, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2 = scattnlay(
        #     np.array([x]), np.array([m]))
        # print("Qabs = ")
        a = 0
        #



crossplane='XZ'
#crossplane='YZ'
#crossplane='XY'

# Options to plot: Eabs, Habs, Pabs, angleEx, angleHy
fields = ['Eabs','Habs']
field_to_plot='Eabs'
#field_to_plot='angleEx'
comment='SiAgSi-absorber-flow'
WL_units='nm'

import os, glob
import string
names = []
for pathAndFilename in glob.iglob(os.path.join(r'.', r'*')):
    title, ext = os.path.splitext(os.path.basename(pathAndFilename))
    if "--nf" in title:
        names.append(title+ext)
        #os.rename(pathAndFilename, "./"+str(i).zfill(3)+"__"+title+ext)
print("Files to plot:")
print(names)
for field_to_plot in fields:
    for filename in sorted(names):
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        rcParams.update({'font.size': 14})
        fig, axs = plt.subplots(1,1)#, sharey=True, sharex=True)
        fig.tight_layout()
        WL=3.75
        fieldplot(fig, axs, WL, filename, comment, WL_units, crossplane, field_to_plot, outline_width=1.5)
        fig.subplots_adjust(hspace=0.3, wspace=-0.1)
        plt.savefig(filename+"-MSTM-"
                        +field_to_plot+".pdf",pad_inches=0.02, bbox_inches='tight')
        plt.draw()
        #    plt.show()
        plt.clf()
        plt.close()


# WL=800
# cst_data_txt="N2__R_200_100__D_310_Sep_10__nf.dat"
# Er, coordZ, coordX = GetCSTField(cst_data_txt,WL)
# print(Er)
