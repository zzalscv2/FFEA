import sys, os

# Only argument should be a .ffea input file
if len(sys.argv) != 2:
	print "Usage: python FFEA_plot_everything [.ffea input file]"
	sys.exit()


# Getting measurement filenames
ffea_script = open(sys.argv[1], "r")
lines = ffea_script.readlines()	
for line in lines:
	if line.strip() == "" or line == "\n" or line.strip() == "<param>" or line.strip() == "</param>" or line.strip() == "<blob>" or line.strip() == "</blob>" or line.strip() == "<system>" or line.strip() == "</system>":
		continue

	line = line.strip().split("=")
	tag = line[0].strip()[1:]
	value = line[1].strip()[:-1]
	
	if tag == "measurement_out_fname":
		measurement_fname_root = value[0:-4]
		measurement_fname_ext = value[-4:]
		print measurement_fname_ext
	if tag == "num_blobs":
		num_blobs = int(value)
		break

world_measurement_fname = measurement_fname_root + "_world" + measurement_fname_ext
blob_measurement_fname = []
for i in range(num_blobs):
	blob_measurement_fname.append(measurement_fname_root + "_blob" + str(i) + measurement_fname_ext)

# Plotting everything using octave
for i in range(num_blobs):
	os.system("octave " + os.path.dirname(sys.argv[0]) + "/FFEA_plot_blob_data.m " + blob_measurement_fname[i] + " " + str(i))

# Plotting world data
os.system("octave " + os.path.dirname(sys.argv[0]) + "/FFEA_plot_world_data.m " + world_measurement_fname + " " + str(num_blobs))

