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
	while(fin.readline().strip() != "@ s0 legend \"Temperature\""):
		pass

	while(True):
		line = fin.readline().strip()
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
plt.ylabel("Temp (K)")
plt.xlabel("Time (ps)")
plt.title("GROMACS Temperature Equilibration")
#plt.axis([time[0] - (time[0] + time[-1]) / 10.0, time[-1], temp.min() * 10, temp.min() * -2])
plt.show()

