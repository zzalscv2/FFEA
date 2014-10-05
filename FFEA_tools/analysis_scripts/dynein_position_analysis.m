arg_list = argv();
data = dlmread(arg_list{1}, "\t");
data = data(2:end,[1 2 4 5]);
frame = data(:, 1);
separation_score = data(:, 2);
position_score = data(:, 3);
orientation_score = data(:, 4);
plot(frame, position_score)
#pause()

separation_bin = zeros(10, 1);
position_bin = zeros(10, 1);
orientation_bin = zeros(10, 1);

total = 0;
for i=1:size(separation_score)
	total += 1;
	for j=0:9
		if separation_score(i) >= j * 0.1 && separation_score(i) < (j + 1) * 0.1
			separation_bin(j + 1) += 1;
		endif
			
		if position_score(i) >= j * 0.1 && position_score(i) < (j + 1) * 0.1
			position_bin(j + 1) += 1;
		endif	

		if orientation_score(i) >= j * 0.1 && orientation_score(i) < (j + 1) * 0.1
			orientation_bin(j + 1) += 1;
		endif
	endfor
endfor

indices = zeros(10, 1);
for i=1:size(indices)
	indices(i) = i/10;
endfor

for i=1:size(position_bin)
	position_bin(i) /= total;
endfor
for i=1:size(separation_bin)
	separation_bin(i) /= total;
endfor
for i=1:size(orientation_bin)
	orientation_bin(i) /= total;
endfor
bar(indices, position_bin)
title("Probability of being in each Position")
xlabel("Normalised Overlap  1 - r_{cm-cm}/r_{cm-cm max}")
ylabel("Probability")
print("position_graph_12.jpg")

bar(indices, orientation_bin)
title("Probability of being in each Orientation")
xlabel("n_{cm1}·n_{cm2}")
ylabel("Probability")
print("orientation_graph_12.jpg")
