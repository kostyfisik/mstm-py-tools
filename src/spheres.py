#!/usr/bin/env python
# -*- coding: UTF-8 -*-
class Spheres:
	radii = [] #nm
	coords = [[], [], []] #nm
	ref_index = []
	def AddSphere(self, radius, coords, index):
		assert (radius > 0),("ERROR! Found R<0")
		self.radii.append(radius)
		for i in range(3):
			self.coords[i].append(coords[i])
		self.ref_index.append(index)
	def Print(self, i):
		return("{:.15E}".format(self.radii[i])+
			       " {:.15E}".format(self.coords[0][i])+
			       " {:.15E}".format(self.coords[1][i])+
			       " {:.15E}".format(self.coords[2][i])+
			       " {:.15E}".format(self.ref_index[i].real)+
			       " {:.15E}".format(self.ref_index[i].imag)
			       )
	def PrintAll(self):
		print_all=""
		for i in range(self.Count()):
			print_all += self.Print(i)+"\n"
		return print_all
	def Count(self):
		return len(self.radii)
# spheres = Spheres()
# spheres.AddSphere(0.12345678901234567890e2, [0, 0.1, 0.2], 2+3j)
# spheres.AddSphere(0.890e2, [10, 10.11, 10.2], 12+13j)
# print(spheres.Count())
# print(spheres.PrintAll())
