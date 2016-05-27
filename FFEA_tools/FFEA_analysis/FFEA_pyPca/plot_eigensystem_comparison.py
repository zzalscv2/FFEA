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
for i in range(num_modes):
	row_labels.append(str(i))
	column_labels.append(str((num_modes - 1) - i))

ax = plt.subplot(111)
ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
ax.set_xticklabels(row_labels)
ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
ax.set_yticklabels(column_labels)
heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
plt.title("Dot Products Between two Eigensystems")
plt.xlabel("EigenSystem 1")
plt.ylabel("Eigensystem 2")
#legend
cbar = plt.colorbar(heatmap)
cbar.set_label('Eigenstate - Eigenstate Dot Products')
plt.savefig(base + ".png")
plt.show()
