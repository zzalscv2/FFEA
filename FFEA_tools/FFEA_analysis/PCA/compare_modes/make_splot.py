import sys, os
import math

if len(sys.argv) != 4:
	sys.exit("Usage: python make_splot.py [LIST OF DOT PRODUCTS] [X DIM] [Y DIM]")

xdim = int(sys.argv[2])
ydim = int(sys.argv[3])
dot_products = open(sys.argv[1], "r")
values = [abs(float(line)) for line in dot_products.readlines()]
dot_products.close()

grid = [values[n:n+xdim] for n in range(0, ydim*xdim, xdim)]

for i in xrange(ydim):
	for j in xrange(xdim):
		v = grid[i][j]
		print i, j, v
		print i, (j + 0.99999), v
	print ""
	for j in xrange(xdim):
		v = grid[i][j]
		print (i + 0.99999), j, v
		print (i + 0.99999), (j + 0.99999), v
	print ""
