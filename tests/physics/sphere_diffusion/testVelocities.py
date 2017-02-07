import sys, os
import FFEA_script, FFEA_trajectory, FFEA_measurement
import numpy as np
from matplotlib import pyplot as plt

script = FFEA_script.FFEA_script(sys.argv[1])
traj = FFEA_trajectory.FFEA_trajectory(script.params.trajectory_out_fname, num_frames_to_read=20000)
meas = FFEA_measurement.FFEA_measurement(script.params.measurement_out_fname, num_frames_to_read=20000)
top = script.load_topology(0)
node = script.load_node(0)
mat = script.load_material(0)
velsum2 = np.array([0.0 for i in range(traj.num_frames)])

i = 0
for f in traj.blob[0][0].frame:
	vs2 = np.array([0.0,0.0,0.0])
	for v in f.vel:
		vs2 += v
	vs2 *= 1.0 / len(top.get_linear_nodes())
	velsum2[i] = np.dot(vs2,vs2)
	i += 1

# m velsum2 == 3kT
# m vel2sum == 3 * num_nodes * kT
vs2 = np.mean(velsum2)
mass = top.calc_mass(mat=mat, node=node, scale = script.blob[0].scale)
print mass
print("mv2 = %e, 3kT = %e\n" % (mass * vs2, 3 * script.params.kT))

EkT = 0.5 * mass * velsum2 * (1.0 / script.params.kT)
EkTavg = np.array([np.mean(EkT[0:i+1]) for i in range(len(EkT))])
plt.plot(meas.global_meas["Time"] * 1e9, EkTavg)
plt.plot(meas.global_meas["Time"] * 1e9, np.array([1.5 for i in range(len(EkT))]))
plt.xlabel("Time (ns)")
plt.ylabel("0.5 m<v^2> / kT")
plt.show()
