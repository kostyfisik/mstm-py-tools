#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import Spheres as sp
import math
import numpy as np
class InputFile:
    ############################################################################
    #Geom
    D = 0
    Sep = 0
    axis = 'x'
    spheres = sp.Spheres()
    ############################################################################
    #Source. The first one is used for ordinary simulaitons
    # spheres.WL = 800
    ang_a = 0.0
    ang_b = 0.0
    ang_pol = 0.0
    # The list of sources can be filled to be iterated over for
    # evaluation of the overlap integral.
    source_WL = []
    #Source angle
    ang_a_kz_x = []
    ang_b_k_z = []
    ang_pol_e_kz = []
    ############################################################################
    #Field
    # isPlotField = False
    isPlotField = True
    sign = ""
    cut_plane = 'xz'
    cut_plane_values={'xy':3, 'yx':3, 'yz':1, 'zy':1, 'zx':2, 'xz':2}
    ############################################################################
    cut_plane = 'xz'
    cut_plane_values={'xy':3, 'yx':3, 'yz':1, 'zy':1, 'zx':2, 'xz':2}
    #plot_scale = 1.0 # Ratio to first sphere
    plot_scale = 0.99995 # Ratio to first sphere for 
    plot_points_per_diameter = 6 # for the first sphere
    nf_plane_position = 0.0
    r_in = 0.0;    r_out = 0.0
    ############################################################################
    def IntegrateOverlapNF(self, in_data, span_value):
        """This will evaluate overlapping integral for field from two sources.
        """
        self.isPlotField = True
        self.plot_scale = 0.99995 # Should be a bit <1, to get points
                                  # on the boundary correct
        out_data = in_data
        R = self.spheres.radii[0]
        plane_pos = np.linspace(-R, R, self.plot_points_per_diameter)
        Es, Hs, coordXs, coordZs  = [], [], [], []
        self.r_in = 0.0;        self.r_out = 0.0;        integral = 0.0
        for pos in plane_pos:
            del Es[:]; del Hs[:]; del coordXs[:]; del coordZs[:]
            for i in range(len(self.source_WL)):
                self.GetSource(i)  #Set WL first
                scale = 2.0*math.pi/float(self.spheres.WL)
                self.nf_plane_position = pos*scale
                self.Run()
                out_data = self.ReadData(out_data, span_value)
                E, H, coordX, coordZ = self.ReadField(self.sign+"--nf.dat",self.spheres.WL)
                Es.append(E);    Hs.append(H)
                coordXs.append(coordX/scale);    coordZs.append(coordZ/scale)
            self.CheckCoords(coordXs)
            self.CheckCoords(coordZs)
            integral += self.IntegrateInPlane(Es, pos, coordXs[0], coordZs[0])
        print ("r_in share = "+str(self.r_in/(self.r_out+self.r_in)))
        return out_data, integral/self.r_in
    ############################################################################
    def IntegrateInPlane(self,Es, coord_y, coordXs, coordZs):
        print(Es[0][0])
        plane_int = 0.0
        for i in range(len(coordXs)):
            dist=math.sqrt(coordXs[i]**2+coord_y**2+coordZs[i]**2)
            if dist < self.spheres.radii[0]*self.plot_scale:
                x = Es[0][i]
                #E0 = [x[0]+1.0j*x[1], x[2]+1.0j*x[3], x[4]+1.0j*x[5]]
                E0 = np.sqrt(np.absolute(x[0]+1.0j*x[1])**2+np.absolute(x[2]+1.0j*x[3])**2+np.absolute(x[4]+1.0j*x[5])**2)
                x = Es[1][i]
                #E1 = [x[0]+1.0j*x[1], x[2]+1.0j*x[3], x[4]+1.0j*x[5]]
                E1 = np.sqrt(np.absolute(x[0]+1.0j*x[1])**2+np.absolute(x[2]+1.0j*x[3])**2+np.absolute(x[4]+1.0j*x[5])**2)
                plane_int += E0*E1
                self.r_in += 1
            else:
                self.r_out +=1
        return plane_int
    ############################################################################
    def CheckCoords(self,coords):
        ref = coords[0]
        for coord in coords:
            if len(coord) != len(ref):
                 raise ValueError("ERROR!!! Cross-sections coordinats do not match by len() !")
            i = 0
            for i in range(len(coord)):
                diff = 1.0 - coord[i]/ref[i]
                if abs(diff) > 1e-2:  # Defined by number of digits in the MSTM output
                    raise ValueError("ERROR!!! Cross-sections coordinats do not match: %f and %f" % (coord[i], ref[i]))
    ############################################################################
    def ReadField(self,data_txt, WL):
        skips = 0
        with open(data_txt, 'r') as data_file:
            for data_line in data_file:
                if len(data_line.split()) > 4: break
                skips += 1
        x=np.transpose(np.loadtxt(data_txt, skiprows=skips))
        E = x[1:7].transpose()
        print(E)
        # E = np.sqrt(np.absolute(x[2]+1.0j*x[3])**2+np.absolute(x[4]+1.0j*x[5])**2+np.absolute(x[6]+1.0j*x[7])**2)
        H = []
        H.append(x[7:13])
        # H = np.sqrt(np.absolute(x[8]+1.0j*x[9])**2+np.absolute(x[10]+1.0j*x[11])**2+np.absolute(x[12]+1.0j*x[13])**2)
        # Coordinates are defined by a cut plane. Here the names are for `xz` case.
        coordX = x[0]
        coordZ = x[1]
        # coordX = np.unique(x[0])
        # coordZ = np.unique(x[1])
        return E, H, coordX, coordZ
    ############################################################################
    def ReadData(self, in_data, span_value):
        out_data = in_data
        isData = False
        WL = self.spheres.WL
        with open(self.sign+'.dat', 'r') as data_file:
            for data_line in data_file:
                if "Qsca" in data_line:
                    if len(out_data) == 0:
                        out_data += "# x WL "+data_line
                    isData = True
                    continue
                if isData == True:
                    # print(data_line[:-1])
                    if 12 == len(data_line.split()):
                        out_data += str(span_value)+" "+str(WL)+" "+data_line
                    else:
                        isData = False
        return out_data
    ############################################################################
    def AddSource(self,WL, ang_a, ang_b, ang_pol):
        self.source_WL.append(WL)
        self.ang_a_kz_x.append(ang_a)
        self.ang_b_k_z.append(ang_b)
        self.ang_pol_e_kz.append(ang_pol)
    ############################################################################
    def ResetSources(self):
        del self.source_WL[:]
        del self.ang_a_kz_x[:]
        del self.ang_b_k_z[:]
        del self.ang_pol_e_kz[:]
    ############################################################################
    def ResetMstmModel(self, WL, R1, R2, n1_name, n2_name, Sep, axis):
        self.Sep = Sep
        self.axis = axis
        self.spheres.Reset()
        self.spheres.WL = WL
        if R1 != 0:
            self.spheres.AddSphere(R1, [0, 0, 0], n1_name)
        if R2 != 0:
            if axis == 'x':
                self.spheres.AddSphere(R2, [R1+R2+Sep, 0, 0], n2_name)
            if axis == 'y':
                self.spheres.AddSphere(R2, [0, R1+R2+Sep, 0], n2_name)
            if axis == 'z+':
                self.spheres.AddSphere(R2, [0, 0, R1+R2+Sep], n2_name)
            if axis == 'z-':
                self.spheres.AddSphere(R2, [0, 0, -(R1+R2+Sep)], n2_name)
            self.SetSign()
    def Run(self):
        self.WriteFile()
        from subprocess import call
        import os
        with  open(os.devnull, 'w') as FNULL:
            call(["mpirun", "-np", "2", "./mstm", "mstm.inp"],
                 stdout=FNULL)
    ############################################################################
    def GetSource(self, i):
        assert (i<len(self.source_WL)),("ERROR! Not enough sources!")
        self.spheres.WL = self.source_WL[i]
        self.ang_a = self.ang_a_kz_x[i]
        self.ang_b = self.ang_b_k_z[i]
        self.ang_pol = self.ang_pol_e_kz[i]
        self.SetSign()
    ############################################################################
    def PrintInput(self):
        self.SetSign()
        scale = 2.0*math.pi/float(self.spheres.WL)
        R1 = self.spheres.radii[0]*scale
        plot_nf = 0
        if self.isPlotField:
            plot_nf = 1
        field_span = R1*self.plot_scale
        step = 2*field_span/float(self.plot_points_per_diameter-1)
        text = """begin_comment
***********************************************************
"""
        text += self.sign + "\nWL = " + str(self.spheres.WL) + str("\n")
        for i in range(self.spheres.Count()):
            text += "R"+str(i)+" = " + str(self.spheres.radii[i]) + str("\n")
        text += "D = "+str(self.D)+"   (Sep = "+str(self.Sep)+""")
***********************************************************
end_comment
number_spheres
"""
        text += str(self.spheres.Count()) + """
sphere_position_file

output_file
"""
        text += self.sign + """.dat
append_output_file
0
run_print_file

write_sphere_data
1
length_scale_factor
"""
# real_chiral_factor
# 0.0d0
# imag_chiral_factor
# 0.0d0
# medium_real_chiral_factor
# 0.d0
# medium_imag_chiral_factor
# 0.d0
        text += "{:.15E}".format(scale) + """
real_ref_index_scale_factor
1.0d0
imag_ref_index_scale_factor
1.0d0
medium_real_ref_index
1.0d0
medium_imag_ref_index
0.d0
target_euler_angles_deg
0.0d0,0.0d0,0.0d0
mie_epsilon
-16
translation_epsilon
1.0d-6
solution_epsilon
1.0d-8
max_number_iterations
5000
store_translation_matrix
0
sm_number_processors
10
near_field_distance
1.0d8
iterations_per_correction
20
incident_or_target_frame
0
normalize_scattering_matrix
1
incident_azimuth_angle_deg
"""
        text += "{:.15E}".format(self.ang_a) + """
incident_polar_angle_deg
"""
        text += "{:.15E}".format(self.ang_b) + """
calculate_scattering_coefficients
1
scattering_coefficient_file

track_iterations
1
calculate_near_field
"""
        text += str(plot_nf) + """
near_field_plane_coord
"""
        text += str(self.cut_plane_values[self.cut_plane]) + """
near_field_plane_position
"""
        text += "{:.15E}".format(self.nf_plane_position) + """
near_field_plane_vertices
"""
        text += "{:.15E} {:.15E} {:.15E} {:.15E}".format(-field_span, -field_span, field_span, field_span) + """
spacial_step_size
"""
        text += "{:.15E}".format(step) + """
polarization_angle_deg
"""
        text += "{:.15E}".format(self.ang_pol) + """
near_field_output_file
"""
        text += self.sign + """--nf.dat
near_field_output_data
2
plane_wave_epsilon
1.0d-4
calculate_t_matrix
1
t_matrix_file
tm_default.dat
t_matrix_convergence_epsilon
1.d-7
sphere_sizes_and_positions
"""
        text += self.spheres.PrintAll()
        return text
    ############################################################################
    def SetSign(self):
        self.D = 0    
        if self.spheres.Count() == 2:
            for i in range(3):
                self.D += (self.spheres.coords[i][0]-self.spheres.coords[i][1])**2
            self.D = math.sqrt(self.D)
            self.Sep = self.D - (self.spheres.radii[0]+self.spheres.radii[1])
            assert (self.Sep>=0), ("ERROR! Found Sep<0")
        self.sign = "N" + str(self.spheres.Count())+"--R-" + str(self.spheres.radii[0])
        if self.spheres.Count() == 2:
            self.sign += "-"+ str(self.spheres.radii[1])
            self.sign += "-D-" + "{:.3g}".format(self.D) +"-"+self.axis+"-Sep-"+"{:.3g}".format(self.Sep)
        self.sign += "-WL-{:04.0f}nm".format(self.spheres.WL)
        if self.isPlotField == True:
            self.sign += "-"+self.cut_plane
    ############################################################################
    def WriteFile(self):
        print(self.PrintInput())
        with open('mstm.inp', 'w') as f:
            f.write(self.PrintInput())
    ############################################################################
    ############################################################################
    ############################################################################
