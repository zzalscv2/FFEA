import sys, os

def check_for(varname):
	if varname in command:
		result = command.replace(varname + '=', '')
		print varname + ' = ' + result
		return result
	return None

if len(sys.argv) != 4:
	sys.exit("Usage: python anal_contact.py [WALRUS FILE] ['ON' DISTANCE] [DUMP RATIO (1 = True, 0 = False)]")

print "Running: anal_contact.py"

path = os.path.dirname(sys.argv[1])
if path == '':
	path = '.'
path += "/"
print "path of walrus file: " + path

command_list = open(sys.argv[1], "r").read()
command_list = ''.join(command_list.split())
command_list = command_list.replace('>', '').split('<')

on_dist = float(sys.argv[2])

dump_ratio = int(sys.argv[3])
if dump_ratio == 1:
	ratiofile = open("__dump_ratio_output__", "w")

vdw_r_eq = None
es_N_x = None
es_N_y = None
es_N_z = None
kappa = None
es_h = None
in_params = False
t_on = 0
t_off = 0
total_steps = 0
num_steps = 0
for command in command_list:
	if 'param' in command:
		in_params = True
	if '/param' in command:
		in_params = False
		vdw_r_eq = float(vdw_r_eq)
		es_N_x = int(es_N_x)
		es_N_y = int(es_N_y)
		es_N_z = int(es_N_z)
		kappa = float(kappa)
		es_h = float(es_h)
		c = es_h * 1.0/kappa
		lx = c * es_N_x
		ly = c * es_N_y
		lz = c * es_N_z
		areaxz = lx * lz
		vol_off = areaxz * (ly - 2 * on_dist)
		print "Box dimensions = (", lx, ly, lz, ")"
		print "Area x-z plane =", areaxz
		print "Distance from surface at which objects are considered 'on' the surface given as", on_dist
		print "Volume of box in which objects are considered 'off' the surface =", vol_off, " (", vol_off * 1e27," nm^3)"
		print "Fraction of total volume that is 'off'", vol_off/(lx * ly * lz)

	if in_params == True:
		if vdw_r_eq == None:
			vdw_r_eq = check_for('vdw_r_eq')
		if es_N_x == None:
			es_N_x = check_for('es_N_x')
		if es_N_y == None:
			es_N_y = check_for('es_N_y')
		if es_N_z == None:
			es_N_z = check_for('es_N_z')
		if kappa == None:
			kappa = check_for('kappa')
		if es_h == None:
			es_h = check_for('es_h')

	if 'measurement' in command:
		measfname = path + command.replace('measurement=', '')
		print "Analysing " + measfname
		measfile = open(measfname, "r")
		for entry in measfile.readlines()[1:]:
			entry_split = entry.split()
			steps = entry_split[0]
			com_y = float(entry_split[4])
			total_steps += 1
			if com_y < on_dist or ly - com_y < on_dist:
				t_on += 1
			else:
				t_off += 1				

			if dump_ratio == 1:
				if t_off > 1:
					ratio = float(t_on)/float(t_off)
					ratiofile.write(str(total_steps) + " " + str(ratio) + "\n")
		measfile.close()
		num_steps = steps

if dump_ratio == 1:
	ratiofile.close()

print "aggregated t_on =", t_on
print "aggregated t_off =", t_off
ratio = float(t_on)/float(t_off)
print "t_on/t_off =", ratio
print "vol_off/(2 * areaxz) * t_on/t_off =", vol_off/(2 * areaxz) * ratio
print "num_simulation_steps (single blob) =", num_steps
print "Done. --> anal_contact.py"
