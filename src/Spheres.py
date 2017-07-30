#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import numpy as np
class Spheres:
    radii = [] #nm
    coords = [[], [], []] #nm
    ref_index = []
    materials = [] #files "BaTiO3-Wemple-o" or "Au-Jhonson" or "Au-Rakic" + ".txt"
    WL = 0
    def GetRefIndex(self,i):
        WLs = np.linspace(self.WL, self.WL, 2)
        fname = self.materials[i]+".txt"
        return self.GetIndex(WLs, fname)[0][1]
    def GetIndex(self,WLs, fname):
        print(fname)
        data = np.loadtxt(fname)
        WL = data[:,0]*1000.0
        indexRe = data[:,1]
        indexIm = np.zeros(data.shape[0])
        if data.shape[1] == 3:  #if material losses are defined in data file
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

    def AddSphere(self, radius, coords, material):
        assert (radius > 0),("ERROR! Found R<0")
        self.radii.append(radius)
        for i in range(3):
            self.coords[i].append(coords[i])
        self.materials.append(material)
    def Print(self, i):
        ref_index = self.GetRefIndex(i)
        return("{:.14E}".format(self.radii[i])+
                   " {:.14E}".format(self.coords[0][i])+
                   " {:.14E}".format(self.coords[1][i])+
                   " {:.14E}".format(self.coords[2][i])+
                   " {:.14E}".format(ref_index.real)+
                   " {:.14E}".format(ref_index.imag)
                   )
    def PrintAll(self):
        print_all=""
        for i in range(self.Count()):
            print_all += self.Print(i)+"\n"
        return print_all
    def Count(self):
        return len(self.radii)
    def Reset(self):
        del self.radii[:]
        for i in range(3):
            del self.coords[i][:]
        del self.ref_index[:]
        del self.materials[:]
        WL = 0
		
# spheres = Spheres()
# spheres.AddSphere(0.12345678901234567890e2, [0, 0.1, 0.2], 2+3j)
# spheres.AddSphere(0.890e2, [10, 10.11, 10.2], 12+13j)
# print(spheres.Count())
# print(spheres.PrintAll())
# spheres.Reset()
# spheres.AddSphere(50.890e2, [510, 510.11, 510.2], 512+513j)
# print(spheres.Count())
# print(spheres.PrintAll())
