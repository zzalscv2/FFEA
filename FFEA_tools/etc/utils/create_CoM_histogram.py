import os, sys, math

if len(sys.argv) != 6:
	sys.exit("Usage: python create_CoM_histogram.py [INPUT MEASUREMENT FILE] [BIN WIDTH] [DIM Y] [kT] [EQUILIBRATION]")

inputfname = sys.argv[1]

bin_width = float(sys.argv[2])

graph_fname = "__create_CoM_histogram__outfile__"
outfile = open(graph_fname, 'w')

outfile.write("# bin width = " + str(bin_width) + "\n")

dim_y = float(sys.argv[3])

kT = float(sys.argv[4])

equilibration = int(sys.argv[5])

outfile.write("# Extracting CoM y values, calculating delta_U and separating into on/off ranges..." + "\n")
infile = open(inputfname, "r")

data_on = []
data_off = []

for line in infile.readlines()[equilibration:]:
	if line.startswith('#'):
		continue
	sl = line.split()
	y = float(sl[4])
	elastic_U = float(sl[2])
	vdw_U = float(sl[12])
	vdw_area = float(sl[10])
	if vdw_area > 0:
		if (dim_y - y) < y:
			y = dim_y - y

		data_on.append([y, elastic_U + vdw_U])
	else:
		data_off.append([y, elastic_U + vdw_U])
outfile.write("# ...done." + "\n")

outfile.write("# Finding max and min values of CoM y..." + "\n")
y_max = max(data_on, key = lambda x: x[0])[0]
y_min = min(data_on, key = lambda x: x[0])[0]

outfile.write("# ...done. y_max = " + str(y_max) + ", y_min = " + str(y_min) + "\n")

num_bins = int(math.ceil((y_max - y_min) / bin_width))
outfile.write("# num bins needed = " + str(num_bins) + "\n")

num_points = len(data_on)
outfile.write("# num_points = " + str(num_points) + "\n")

outfile.write("# Binning..." + "\n")
hist = [[0, 0.0] for i in range(num_bins)]

for entry in data_on:
	index = int(math.floor((entry[0] - y_min)/bin_width))
	(hist[index])[0] += 1
	(hist[index])[1] += entry[1]

for i in range(len(hist)):
	if hist[i][0] != 0:
		hist[i][1] /= hist[i][0]
	hist[i][0] = float(hist[i][0])/num_points
outfile.write("# ...done." + "\n")

k = 1.3806503e-23
T = kT/k

outfile.write("# kT = " + str(kT) + "\n")
outfile.write("# k = " + str(k) + "\n")
outfile.write("# T = " + str(T) + "\n")

outfile.write("# CoM_y | prob | average dU | dF | dS" + "\n")

for i in range(num_bins):
	y = i * bin_width + y_min
	p = (hist[i])[0]
	dU = (hist[i])[1]
	if p > 0.0:
		dF = -kT * math.log(p)
		dS = (dU - dF)/T
	else:
		dF = float("NaN")
		dS = float("NaN")
	outfile.write(str(y) + " " + str(p) + " " + str(dU) + " " + str(dF) + " " + str(dS) + "\n")

outfile.close()

gnuplotscriptfname = "__create_CoM_histogram_gnuplotscript__"
print "gnuplot script will be " + gnuplotscriptfname
gnuplotscriptfile = open(gnuplotscriptfname, "w")
gnuplotscriptfile.write("set term postscript eps enhanced color font 'Helvetica,8'\n")
gnuplotscriptfile.write('set output "create_CoM_histogram_all.eps"\n')
gnuplotscriptfile.write("set multiplot layout 2,2 rowsfirst\n")

gnuplotscriptfile.write("unset key\n")
gnuplotscriptfile.write('set xlabel "Distance from surface"\n')
gnuplotscriptfile.write('set ylabel "Probability"\n')
gnuplotscriptfile.write("plot './" + graph_fname + "' u 1:2 w lp\n")

gnuplotscriptfile.write("unset key\n")
gnuplotscriptfile.write('set xlabel "Distance from surface"\n')
gnuplotscriptfile.write('set ylabel "dU = U_vdw + U_elastic"\n')
gnuplotscriptfile.write("plot './" + graph_fname + "' u 1:3 w lp\n")

gnuplotscriptfile.write("unset key\n")
gnuplotscriptfile.write('set xlabel "Distance from surface"\n')
gnuplotscriptfile.write('set ylabel "dF = -kT ln(p)"\n')
gnuplotscriptfile.write("plot './" + graph_fname + "' u 1:4 w lp\n")

gnuplotscriptfile.write("unset key\n")
gnuplotscriptfile.write('set xlabel "Distance from surface"\n')
gnuplotscriptfile.write('set ylabel "dS = (dU - dF)/T"\n')
gnuplotscriptfile.write("plot './" + graph_fname + "' u 1:5 w lp\n")

gnuplotscriptfile.write("unset multiplot\n")
gnuplotscriptfile.close()

os.system("gnuplot " + gnuplotscriptfname + "\n")
