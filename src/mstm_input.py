#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import Spheres as sp
import math
class InputFile:
	spheres = sp.Spheres()
	WL = 800
	D = 0
	Sep = 0
	sign = ""
	isPlotField = False
	# cut_plane = 'xy'
	cut_plane = 'xz'
	cut_plane_values={'xy':3, 'yx':3, 'yz':1, 'zy':1, 'zx':2, 'xz':2}
	plot_scale = 1
	points = 501
	def PrintInput(self):
		self.SetSign()
		scale = 2.0*math.pi/float(self.WL)
		R1 = self.spheres.radii[0]*scale
		plot_nf = 0
		if self.isPlotField: plot_nf = 1
		span = round(2.0*R1*self.plot_scale*self.points)/self.points
		step = 2*span/float(self.points)
		
		text = """begin_comment
***********************************************************
"""
		text += self.sign + "\nWL = " + str(self.WL) + str("\n")
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
0.d0
incident_polar_angle_deg
0.0d0
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
0.0d0
near_field_plane_vertices
"""
		text += "{:f} {:f} {:f} {:f}".format(-span, -span, span, span) + """
spacial_step_size
"""
		text += "{:f}".format(step) + """
polarization_angle_deg
0.d0
near_field_output_file
"""
		text += self.sign + """__nf.dat
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
	def SetSign(self):
		self.D = 0		
		if self.spheres.Count() == 2:
			for i in range(3):
				self.D += (self.spheres.coords[i][0]-self.spheres.coords[i][1])**2
			self.D = math.sqrt(self.D)
			self.Sep = self.D - (self.spheres.radii[0]+self.spheres.radii[1])
			assert (self.Sep>=0), ("ERROR! Found Sep<0")
		self.sign = "N" + str(self.spheres.Count())+"__R_" + str(self.spheres.radii[0])
		if self.spheres.Count() == 2:
			self.sign += "_"+ str(self.spheres.radii[1])
			self.sign += "__D_" + "{:.3g}".format(self.D) +"_Sep_"+"{:.3g}".format(self.Sep)
		self.sign += "__WL_{:04.0f}nm".format(self.WL)
                if self.isPlotField == True:
                    self.sign += "__"+self.cut_plane

	def WriteFile(self):
		with open('mstm.inp', 'w') as f:
			f.write(self.PrintInput())
			

# R1 = 200
# n1 = 2+0.1j
# n2 = 3+2j
# R2 = 100
# Sep = 10

# mstm_input = InputFile()
# mstm_input.spheres.AddSphere(R1, [0, 0, 0], n1)
# mstm_input.spheres.AddSphere(R2, [R1+R2+Sep, 0, 0], n2)
# mstm_input.cut_plane = 'xy'
# mstm_input.plot_scale = (R1+3*Sep+2*R2)/R1
# mstm_input.points = 150
# #mstm_input.isPlotField = False
# mstm_input.isPlotField = True

# print(mstm_input.PrintInput())
# mstm_input.WriteFile()
