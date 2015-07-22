import sys
from math import radians
from math import cos as c
from math import sin as s

if len(sys.argv) != 4:
	sys.exit("Usage: python Generate_rotation_matrix.py [x angle] [y angle] [z angle] in degrees!")

print "Rotation matrix produced: R = R(z)R(y)R(x) i.e. x axis, then y, then z\n"

# Get args
x = radians(float(sys.argv[1]))
y = radians(float(sys.argv[2]))
z = radians(float(sys.argv[3]))

# Build matrix
M = [0.0 for i in range(9)]
M[0] = c(y) * c(z)
M[1] = s(x) * s(y) * c(z) - c(x) * s(z)
M[2] = c(x) * s(y) * c(z) + s(x) * s(z)
M[3] = c(y) * s(z)
M[4] = s(x) * s(y) * s(z) + c(x) * c(z)
M[5] = c(x) * s(y) * s(z) - s(x) * c(z)
M[6] = -1 * s(y)
M[7] = s(x) * c(y)
M[8] = c(x) * c(y)

# Print it
print M[0], M[1], M[2]
print M[3], M[4], M[5]
print M[6], M[7], M[8]



