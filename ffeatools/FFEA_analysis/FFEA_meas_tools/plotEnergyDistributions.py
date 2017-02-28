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

import sys, os
import FFEA_script, FFEA_measurement
import argparse as _argparse
import __builtin__

try:
    from matplotlib import pyplot as plt
    import numpy as np
except(ImportError):
    raise ImportError("ImportError recieved. Without the numpy and matplotlib libraries, this script cannot run. Please install these, or manually use gnuplot :)")

# Set up argparse
parser = _argparse.ArgumentParser(description="Plot the probability distribution of energy for your system")
parser.add_argument("script", action="store", help="Input script file (.ffea).")

def plotEnergyDistributions(script):
	
	meas = script.load_measurement() 
	kT = script.params.kT

	#
	# First, global graph (strain and kinetic)
	#
	cmeas = meas.global_meas
	num_steps = len(cmeas["Time"])

	# Get energy data (kbt)
	Estrain = cmeas["StrainEnergy"] * (1.0 / kT)
	Estrain = Estrain[500:]
	Estrainmean = np.mean(Estrain)
	Estrainstddev = np.std(Estrain)
	#for i in reversed(range(len(Estrain))):
	#	if Estrain[i] > Estrainmean + 3 * Estrainstddev or Estrain[i] < Estrainmean - 3 * Estrainstddev:
	#		Estrain = np.delete(Estrain, i)

	if script.params.num_blobs >= 1:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.hist(Estrain, 50, normed = 1, facecolor='g', alpha=0.75, label = "Simulation")
	plt.show()

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    # Get args and build objects
    script = FFEA_script.FFEA_script(args.script)
    if not os.path.exists(script.params.measurement_out_fname):
        raise IOError("Error. Measurement file not found. Please supply a script with a completed measurement file.")
    plotEnergyDistributions(script)
