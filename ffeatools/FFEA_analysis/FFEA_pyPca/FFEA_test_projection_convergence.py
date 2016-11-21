import sys, os
import numpy as np
from matplotlib import pyplot as plt

if len(sys.argv) != 2:
	sys.exit("Usage: python FFEA_get_PCA_projections.py [INPUT .proj file]")

# Get args
pfname = sys.argv[1]
base, ext = os.path.splitext(os.path.abspath(pfname))

# Load projections in
proj = []
fin = open(pfname, "r")
for line in fin:
	proj.append([float(i) for i in line.split()])
fin.close()

for p in proj:
	p = np.array(p)
	pm = np.mean(p)
	psd = np.std(p)
	print pm, psd
	corr = np.correlate(p,p, mode='full')
	corr = corr[corr.size/2:] * 1.0 / np.power(psd,2)
	plt.plot(corr)
	plt.show()
