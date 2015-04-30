import sys, os
import math

if len(sys.argv) != 3:
	sys.exit("Usage: python dot.py [EIGEN VECTOR 1] [EIGEN VECTOR 2]")

eigen1 = open(sys.argv[1], "r")
eigen2 = open(sys.argv[2], "r")

N1 = int((eigen1.readline().split())[2])
N2 = int((eigen2.readline().split())[2])

if N1 != N2:
	sys.exit("Vectors are not of equal size (N1 = " + str(N1) + " and N2 = " + str(N2) + " )")

sum = 0
for i in xrange(N1):
	e1 = float(eigen1.readline())
	e2 = float(eigen2.readline())
	sum += e1 * e2

print sum
