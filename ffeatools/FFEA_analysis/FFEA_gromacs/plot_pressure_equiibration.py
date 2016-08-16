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
temp = []
with open(infname) as fin:
	line = fin.readline().strip()
	while(line != "@ s0 legend \"Pressure\""):
		print line
		line = fin.readline().strip()
		pass

	while(True):
		line = fin.readline().strip()
		print line
		if line == "":
			break
		else:
			sline = line.split()
			time.append(float(sline[0]))
			temp.append(float(sline[1]))

time = np.array(time)
temp = np.array(temp)

# Plot a pretty graph
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
line, = ax.plot(time, temp, color='blue', lw=2)
#ax.set_yscale('symlog')
plt.ylabel("Pressure (bar)")
plt.xlabel("Time (ps)")
plt.title("GROMACS Pressure Equilibration")
#plt.axis([time[0] - (time[0] + time[-1]) / 10.0, time[-1], temp.min() * 10, temp.min() * -2])
plt.show()

