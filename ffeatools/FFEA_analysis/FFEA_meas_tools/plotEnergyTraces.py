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
import FFEA_script, FFEA_measurement
import argparse as _argparse
import __builtin__

try:
    from matplotlib import pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    import numpy as np
except(ImportError):
    raise ImportError("ImportError recieved. Without the numpy and matplotlib libraries, this script cannot run. Please install these, or manually use gnuplot :)")

# Set up argparse
parser = _argparse.ArgumentParser(description="Test KE of blob against classical equipartition theorem")
parser.add_argument("script", action="store", help="Input script file (.ffea).")

def plotEnergyTraces(script):

    meas = FFEA_measurement.FFEA_measurement(script.params.measurement_out_fname)
    top = [script.load_topology(i) for i in range(script.params.num_blobs)]
    
    # We need to plot a global measurement graph, and a graph for every blob
    total_num_nodes = 0
    total_num_mass_nodes = 0
    
    for i in range(script.params.num_blobs):
        total_num_nodes += len(top[i].get_linear_nodes())
        if script.blob[i].solver != "CG_nomass":
            total_num_mass_nodes += len(top[i].get_linear_nodes())
    
    kT = script.params.kT
    
    #
    # First, global graph (strain and kinetic)
    #
    cmeas = meas.global_meas
    num_steps = len(cmeas["Time"])

    # Get x axis data (time)
    x = cmeas["Time"] * 1e9    # ns

    if script.params.num_blobs >= 1:
	    fig, ax = plt.subplots()
	    
	    # And y axis data (all energies)
	    ys = []
	    yk = []
	    for i in range(num_steps):
		ys.append(np.mean(cmeas["StrainEnergy"][0:i + 1]) / kT)
		if cmeas["KineticEnergy"] != None:
		    yk.append(np.mean(cmeas["KineticEnergy"][0:i + 1]) / kT)
	    
	    print "\nGlobal System:\n"
	    
	    ysEXP = [(3 * total_num_nodes - 6) / 2.0 for i in range(num_steps)]
	    ysh, = plt.loglog(x, ys, label='ysh')
	    ysEXPh, = plt.loglog(x, ysEXP, "-", label='ysEXPh')
	    
	    yserr = (np.fabs(ys[-1] - ysEXP[-1]) / ysEXP[-1]) * 100.0
	    print "\tTheoretical strain energy = %f; Simulation strain energy = %f; Error is %f%%" % (ysEXP[-1], ys[-1], yserr) 
	    if cmeas["KineticEnergy"] != None:
		ykEXP = [(3 * total_num_mass_nodes) / 2.0 for i in range(num_steps)]
		ykh, = plt.loglog(x, yk, label='ykh')
		ykEXPh, = plt.loglog(x, ykEXP, "-", label='ykEXPh')
	    
		ykerr = (np.fabs(yk[-1] - ykEXP[-1]) / ykEXP[-1]) * 100.0
		print "\tTheoretical kinetic energy = %f; Simulation kinetic energy = %f; Error is %f%%" % (ykEXP[-1], yk[-1], ykerr)
	    
	    # Axis stuff
	  #  ax.set_xlim(xmin=1)

	    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
	    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
	    ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
	    ax.set_xlabel("Time (ns)")
	    ax.set_ylabel("Energy (kT)")
	    ax.set_title("Global Energy - Running Average")
	    
	    # Put details on the graph
	    if cmeas["KineticEnergy"] != None:
		plt.legend([ysh, ysEXPh, ykh, ykEXPh], ['Global Strain Energy - Sim', 'Global Strain Energy - Theory', 'Global Kinetic Energy - Sim', 'Global Kinetic Energy - Theory'], loc = 4)
	    else:
		plt.legend([ysh, ysEXPh], ['Global Strain Energy - Sim', 'Global Strain Energy - Theory'], loc = 4)
	    plt.show()
    
    #
    # Now, individual blobs (strain and kinetic), if we can
    #
    if meas.blob_meas == []:
        sys.exit()

    # Get figure size
    blobs_to_plot = 0
    for i in range(script.params.num_blobs):
	
	if script.blob[i].conformation[0].motion_state == "DYNAMIC":
		blobs_to_plot += 1

    # Too many??
    if blobs_to_plot > 10:
	sys.exit("Too many blobs to print out like this!")

    plot_index = 0
    if blobs_to_plot == 2:
	fig = plt.figure(1, figsize=(10,10))
	subplot_x = 2
	subplot_y = 1
    else:
	fig = plt.figure(1, figsize=(20,20))
	subplot_x = np.ceil(np.sqrt(blobs_to_plot))
    	subplot_y = np.ceil(np.sqrt(blobs_to_plot))
    for i in range(script.params.num_blobs):
	
	if script.blob[i].conformation[0].motion_state == "STATIC":
		print "Skipping STATIC blob " + str(i)
		continue
        
	plot_index += 1
        ax = fig.add_subplot(subplot_x, subplot_y, plot_index)
    
        cmeas = meas.blob_meas[i]
        num_nodes = len(top[i].get_linear_nodes())
    
        # Already got the x axis
    
        # And y axis data (all energies)
        ys = []
        yk = []
        for j in range(num_steps):
            ys.append(np.mean(cmeas["StrainEnergy"][0:j + 1]) / kT)
            if cmeas["KineticEnergy"] != None:
                yk.append(np.mean(cmeas["KineticEnergy"][0:j + 1]) / kT)
    
        print "\nBlob %d:\n" % (i)
        ysEXP = [(3 * num_nodes - 6) / 2.0 for j in range(num_steps)]
        ysh, = ax.loglog(x, ys, label='ysh')
        ysEXPh, = ax.loglog(x, ysEXP, "-", label='ysEXPh')
    
        yserr = (np.fabs(ys[-1] - ysEXP[-1]) / ysEXP[-1]) * 100.0
        print "\tTheoretical strain energy = %f; Simulation strain energy = %f; Error is %f%%" % (ysEXP[-1], ys[-1], yserr) 
        if cmeas["KineticEnergy"] != None:
            ykEXP = [(3 * num_nodes) / 2.0 for j in range(num_steps)]
            ykh, = ax.loglog(x, yk, label='ykh')
            ykEXPh, = ax.loglog(x, ykEXP, "-", label='ykEXPh')
    
            ykerr = (np.fabs(yk[-1] - ykEXP[-1]) / ykEXP[-1]) * 100.0
            print "\tTheoretical kinetic energy = %f; Simulation kinetic energy = %f; Error is %f%%" % (ykEXP[-1], yk[-1], ykerr)

	ax.set_xlim(xmin=1)
  	ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
	ax.yaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
	ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Energy (kT)")
        ax.set_title("Blob %d Energy - Running Average" % (i))
    
        # Put details on the graph
        if cmeas["KineticEnergy"] != None:
            ax.legend([ysh, ysEXPh, ykh, ykEXPh], ['Blob %d Strain Energy - Sim' % (i), 'Blob %d Strain Energy - Theory' % (i), 'Blob %d Kinetic Energy - Sim' % (i), 'Blob %d Kinetic Energy - Theory' % (i)], loc = 4)
        else:
            ax.legend([ysh, ysEXPh], ['Blob %d Strain Energy - Sim' % (i), 'Blob %d Strain Energy - Theory' % (i)], loc = 4)
    plt.show()

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    # Get args and build objects
    script = FFEA_script.FFEA_script(args.script)
    if not os.path.exists(script.params.measurement_out_fname):
        raise IOError("Error. Measurement file not found. Please supply a script with a completed measurement file.")
    plotEnergyTraces(script)
