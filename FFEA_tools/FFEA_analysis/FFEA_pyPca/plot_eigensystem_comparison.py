import sys
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

if len(sys.argv) != 2:
	sys.exit("Usage: python plot_matrix.py [INPUT modes matrix file (.modes)]")

# Get args
modefiles = sys.argv[1]

# Open file and read
fin = open(modefiles, "r")
for i in range(4):
	fin.readline()

data = []
for i in range(5):
	data_line = []
	sline = fin.readline().split()
	for bit in sline:
		data_line.append(abs(float(bit)))
	data.append(data_line)

data = np.array(data)
print data

# Build figure
column_labels = list('01234')
row_labels = list('01234')
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
plt.show()
