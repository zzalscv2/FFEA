import sys, os
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 3:
	sys.exit("Usage: python " + os.path.basename(sys.argv[0]) + " [INPUT .xvg fname] [OUTPUT .jpg fname]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]

# Get data
time = []
energy = []
with open(infname) as fin:
	while(fin.readline().strip() != "@ s0 legend \"Potential\""):
		pass

	while(True):
		line = fin.readline().strip()
		if line == "":
			break
		else:
			sline = line.split()
			time.append(float(sline[0]))
			energy.append(float(sline[1]))

time = np.array(time)
energy = np.array(energy)

# Plot a pretty graph
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
line, = ax.plot(time, energy, color='blue', lw=2)
ax.set_yscale('symlog')
plt.ylabel("Energy (kJ / mol)")
plt.xlabel("Time (ps)")
plt.title("GROMACS Energy Minimisation")
plt.axis([time[0] - (time[0] + time[-1]) / 10.0, time[-1], energy.min() * 10, energy.min() * -2])
plt.show()

