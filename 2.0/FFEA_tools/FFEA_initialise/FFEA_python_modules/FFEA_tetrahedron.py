import sys, os
import numpy as np
import FFEA_triangle

class FFEA_tetrahedron:
	
	def __init__(self, n1, n2, n3, n4):
		
		self.n = [n1,n2,n3,n4]
		self.face = []
		self.face.append(FFEA_triangle.FFEA_triangle(self.n[0] - 1, self.n[1] - 1, self.n[2] - 1))
		self.face.append(FFEA_triangle.FFEA_triangle(self.n[0] - 1, self.n[2] - 1, self.n[3] - 1))
		self.face.append(FFEA_triangle.FFEA_triangle(self.n[0] - 1, self.n[3] - 1, self.n[1] - 1))
		self.face.append(FFEA_triangle.FFEA_triangle(self.n[1] - 1, self.n[3] - 1, self.n[2] - 1))
		
	def calc_volume(self, node_list):
		
		edge = [node_list[self.n[1]] - node_list[self.n[0]], node_list[self.n[2]] - node_list[self.n[0]], node_list[self.n[3]] - node_list[self.n[0]]]
		self.volume = (1.0 / 6.0) * np.linalg.norm(np.dot(np.cross(edge[0], edge[1]), edge[2]))
		
		return self.volume
