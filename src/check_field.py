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
###############################################################################

def GetCSTField(cst_data_txt, WL):
    x=np.transpose(np.loadtxt(cst_data_txt, skiprows=4))
    E = np.sum(np.absolute(x[2]+1.0j*x[3])**2+np.absolute(x[4]+1.0j*x[5])**2+np.absolute(x[6]+1.0j*x[6])**2)
    coordX = np.unique(x[0])/WL
    coordZ = np.unique(x[1])/WL
    print(x[2])
    return E, coordX, coordZ

WL=800
cst_data_txt="N2__R_200_100__D_310_Sep_10__nf.dat"
Er, coordZ, coordX = GetCSTField(cst_data_txt,WL)
print(Er)
