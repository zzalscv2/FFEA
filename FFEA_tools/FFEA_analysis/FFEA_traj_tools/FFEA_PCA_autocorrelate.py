def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]

import sys, os
import numpy as np
from matplotlib import pyplot as plt

if len(sys.argv) != 3:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT PCZ fname] [Projection Number]")

# Get args
pczfname = sys.argv[1]
base, ext = os.path.splitext(sys.argv[1])
projno = int(sys.argv[2])
projfname = base + "_proj%d.out" % (projno)
 
# Get a projection file
os.system("pyPczdump -i %s -p %d -o %s" % (pczfname, projno, projfname))

# Read projections
proj = []
fin = open(projfname, "r")
for line in fin:
	proj.append(float(line))
fin.close()

# Autocorrelate and make graph
acorr = autocorr(proj)
plt.plot(acorr)
plt.show()
