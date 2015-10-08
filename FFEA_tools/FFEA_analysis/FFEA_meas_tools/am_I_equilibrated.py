import sys, os
import FFEA_script, FFEA_measurement, FFEA_topology, FFEA_pin

if len(sys.argv) != 2 and len(sys.argv) != 5:
	sys.exit("Usage: python am_I_equilibrated.py [INPUT .ffea file] OPTIONAL: [num_frames_to_check] [step_length] [kT]")

# Get script
infile = sys.argv[1]

num_frames = float("inf")
if len(sys.argv) == 5:
	num_frames = int(sys.argv[2])

# Open script and load relevant stuff
script = FFEA_script.FFEA_script(infile)
meas = script.load_measurement(num_frames)
topology = [[script.load_topology(i, j) for j in range(script.params.num_conformations[i])] for i in range(script.params.num_blobs)]
pinned = [[script.load_pin(i, j) for j in range(script.params.num_conformations[i])] for i in range(script.params.num_blobs)]

# Get linear nodes - pinned nodes
linear_nodes = [[top.get_num_linear_nodes() for top in blob] for blob in topology]
pinned_nodes = []
for i in range(len(pinned)):
	pinned_set = []
	for j in range(len(pinned[i])):
		pinned_set.append(pinned[i][j].get_num_linear_pinned_nodes(topology[i][j]))

	pinned_nodes.append(pinned_set)

for i in range(len(linear_nodes)):
	for j in range(len(linear_nodes[i])):
		linear_nodes[i][j] -= pinned_nodes[i][j]

if len(sys.argv) == 2:
	meas.plot_energies()
else:
	meas.plot_energies(kT = float(sys.argv[4]), step_length = float(sys.argv[3]), num_nodes = linear_nodes)
