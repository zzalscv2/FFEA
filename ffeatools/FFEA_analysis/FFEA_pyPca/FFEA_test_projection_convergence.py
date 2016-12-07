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
