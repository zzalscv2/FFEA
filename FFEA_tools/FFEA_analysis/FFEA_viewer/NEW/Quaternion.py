import math

class Quaternion:
	def __init__(self):
		self.w = 1
		self.x = 0
		self.y = 0
		self.z = 0

	def rotate(self, angle, axis_x, axis_y, axis_z):
		# Firstly construct the rotation quaternion
		angle *= .5;
		rot_w  = math.cos(angle)
		rot_x = axis_x * math.sin(angle)
		rot_y = axis_y * math.sin(angle)
		rot_z = axis_z * math.sin(angle)

		self.premultiply(rot_w, rot_x, rot_y, rot_z)

		mag = math.sqrt(self.w*self.w +self.x*self.x + self.y*self.y + self.z*self.z)
		self.w /= mag
		self.x /= mag
		self.y /= mag
		self.z /= mag

	# premultiply this quaternion with the given quaternion
	def premultiply(self, pw, px, py, pz):
		self.w = pw * self.w - px * self.x - py * self.y - pz * self.z
		self.x = pw * self.x + px * self.w + py * self.z - pz * self.y
		self.y = pw * self.y - px * self.z + py * self.w + pz * self.x
		self.z = pw * self.z + px * self.y - py * self.x + pz * self.w

	def construct_matrix(self):
		m =	[
			1 - 2*self.y*self.y - 2*self.z*self.z,
			2*self.x*self.y + 2*self.z*self.w,
			2*self.x*self.z - 2*self.y*self.w,
			0,
			2*self.x*self.y - 2*self.z*self.w,
			1 - 2*self.x*self.x - 2*self.z*self.z,
			2*self.y*self.z + 2*self.x*self.w,
			0,
			2*self.x*self.z + 2*self.y*self.w,
			2*self.y*self.z - 2*self.x*self.w,
			1 - 2*self.x*self.x - 2*self.y*self.y,
			0,
			0,
			0,
			0,
			1
			]
		return m

