from os import path
from time import sleep
import FFEA_frame from FFEA_frame

class FFEA_node(FFEA_frame):

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def self.reset()
