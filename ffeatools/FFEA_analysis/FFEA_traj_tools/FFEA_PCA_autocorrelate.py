# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

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
