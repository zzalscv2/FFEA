# Only argument should be a blob measurement file

if( nargin != 2)
	printf("Usage: octave FFEA_plot_blob_data.m [.ffea FILE] [blob_number]")
endif

arg_list = argv();

M = dlmread(arg_list{1});
blob_number = arg_list{2};

# Sparate data into columns
frame = M(3:end,1);
kinetic_energy = M(3:end,2);
potential_energy = M(3:end,3);
rmsd = M(3:end,10);
vdw_xz_area = M(3:end,11);
vdw_xz_force = M(3:end,12);
vdw_xz_energy = M(3:end,13);

# Plot each graph
plot(frame, kinetic_energy);
ylim([0 max(ylim)])
xlabel("Step");
ylabel("Kinetic Energy (J)");
title("Kinetic Energy Trace");
graph_fname = strcat("blob", blob_number, "_kinetic_energy.jpg");
print(graph_fname);

plot(frame, potential_energy);
ylim([0 max(ylim)])
xlabel("Step");
ylabel("Potential Energy (J)");
title("Potential Energy Trace");
graph_fname = strcat("blob", blob_number, "_potential_energy.jpg");
print(graph_fname);

plot(frame, rmsd);
xlabel("Step");
ylabel("rmsd (m^2)");
title("RMSD Trace");
graph_fname = strcat("blob", blob_number, "_rmsd.jpg");
print(graph_fname);

plot(frame, vdw_xz_area);
xlabel("Step");
ylabel("vdw_xz_area (m^2)");
title("Area engaging in vdw interations with xz surface");
graph_fname = strcat("blob", blob_number, "_vdw_xz_area.jpg");
print(graph_fname);

plot(frame, vdw_xz_force);
xlabel("Step");
ylabel("vdw_xz_force (N)");
title("Vdw force with xz surface");
graph_fname = strcat("blob", blob_number, "_vdw_xz_force.jpg");
print(graph_fname);

plot(frame, vdw_xz_energy);
xlabel("Step");
ylabel("vdw_xz_energy (J)");
title("Vdw energyy with xz surface");
graph_fname = strcat("blob", blob_number, "_vdw_xz_energy.jpg");
print(graph_fname);
