import sys, os
import numpy as np

class FFEA_triangle:
	
	def __init__(self, n1, n2, n3):
		
		self.n = [n1,n2,n3]
		
	def calc_area(self, node_list):
		
		edge = np.array([node_list[self.n[1]] - node_list[self.n[0]], node_list[self.n[2]] - node_list[self.n[0]]])
		self.area = 0.5 * np.linalg.norm(np.cross(edge[0], edge[1]))

		return self.area
		
		
