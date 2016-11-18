import sys, os
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

if len(sys.argv) != 2:
	sys.exit("Usage: python plot_matrix.py [INPUT modes matrix file (.modes)]")

# Get args
modefiles = sys.argv[1]
base, ext = os.path.splitext(modefiles)

# Open file and read
fin = open(modefiles, "r")
for i in range(4):
	fin.readline()

data = []
data2 = []
#for i in range(5):

num_modes = 0

while True:
	line = fin.readline()
	if line.strip() == "":
		break
	
	num_modes += 1
	data_line = []
	sline = line.split()
	for bit in sline:
		data_line.append(abs(float(bit)))
		
	
	if num_modes > 1:
		if num_modes_check != len(data_line):
			sys.exit("Error! Expected a square matrix. Make sure you have equal number of eigenvectors in each set.")
	else:
		num_modes_check = len(data_line)
			
	data.append(data_line)
	data2.append(data_line)

for i in range(num_modes):
	data[i] = data2[(num_modes - 1) - i]

data = np.array(data)

# Build figure
#column_labels = list('43210')
#row_labels = list('01234')
column_labels = []
row_labels = []
#xticks = []	# cols
#yticks = []	# rows

# We only want 11 labels if num_modes > 20!

if num_modes <= 20:
	for i in range(num_modes):
		row_labels.append(str(i))
		#column_labels.append(str((num_modes - 1) - i))
		column_labels.append(str(i))
else:
	count = 0
	for i in range(num_modes + 1):
		
		if i % int(num_modes / 10.0) == 0:
			row_labels.append(int((count / 10.0) * num_modes))
			#column_labels.append(int((1 - (count / 10.0)) * num_modes))
			column_labels.append(int((count / 10.0) * num_modes))
			count += 1
		else:
			pass

# Get figure properties
fig, ax = plt.subplots(figsize=(13,10))

# Plot data
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
heatmap.set_clim([0, 1])

# And 'plot' colorbar
cbar = plt.colorbar(heatmap, ticks=[i / 10.0 for i in range(11)])
	
# All tick
if num_modes <= 20:
	ax.set_xticks(np.arange(data.shape[0]) + 0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[1], 0, -1) - 0.5, minor=False)
else:
	ax.set_xticks(np.array([int((i / 10.0) * num_modes) for i in range(11)]), minor=False)
	ax.set_yticks(np.array([int((1 - (i / 10.0)) * num_modes) for i in range(11)]), minor=False)

# Labels
ax.set_yticklabels(row_labels, fontsize=18)
ax.set_xticklabels(column_labels, fontsize=18)

# And now titles and stuff
ax.set_xlim(0, num_modes)
ax.set_ylim(0, num_modes)

ax.set_title("Eigenvector Dot Product Array", fontsize=24)
ax.set_xlabel("Eigensystem 1 Modes", fontsize=18)
ax.set_ylabel("Eigensystem 2 Modes", fontsize=18)

cbar.ax.set_yticklabels(["%2.1f" % (i / 10.0) for i in range(11)], fontsize=18)
cbar.set_label('Normalised Dot Product Values', fontsize=18)

plt.savefig(base + ".png")
plt.show()

#fig, ax = plt.subplots()

#data = np.clip(randn(250, 250), -1, 1)

#cax = ax.imshow(data, interpolation='nearest', cmap=cm.coolwarm)
#ax.set_title('Gaussian noise with vertical colorbar')

# Add colorbar, make sure to specify tick locations to match desired ticklabels
#cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar
