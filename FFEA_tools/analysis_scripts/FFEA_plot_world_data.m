# Only argument should be a blob measurement file

if( nargin != 2)
	printf("Usage: octave FFEA_plot_blob_data.m [.ffea FILE] [num_blobs]")
endif

arg_list = argv();

M = dlmread(arg_list{1});
num_blobs = str2num(arg_list{2});

# Sparate data into columns

frame = M(3:end,1);
t = 1;
for i = 0:num_blobs - 2
	for j = (i + 1):num_blobs-1
		t += 1;

		# Plot area graph
		plot(frame, M(3:end,t));
		xlabel("Step");
		ylabel(strcat("Vdw Area (m^2)"));
		title(strcat("Interacting Vdw Area between Blob", int2str(i), int2str(j)));
		graph_fname = strcat("blob", int2str(i), int2str(j), "_vdw_area.jpg");
		print(graph_fname);

		t += 1;

		# Plot force graph
		plot(frame, M(3:end,t));
		xlabel("Step");
		ylabel(strcat("Vdw Force (N)"));
		title(strcat("Vdw Force between Blob", int2str(i), int2str(j)));
		graph_fname = strcat("blob", int2str(i), int2str(j), "_vdw_force.jpg");
		print(graph_fname);

		t += 1;

		# Plot energy graph
		plot(frame, M(3:end,t));
		xlabel("Step");
		ylabel(strcat("Vdw Energy (J)"));
		title(strcat("Vdw Energy between Blob", int2str(i), int2str(j)));
		graph_fname = strcat("blob", int2str(i), int2str(j), "_vdw_energy.jpg");
		print(graph_fname);
	endfor
endfor

