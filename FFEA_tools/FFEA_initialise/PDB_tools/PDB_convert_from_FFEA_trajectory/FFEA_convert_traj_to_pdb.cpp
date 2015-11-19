#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;
int main(int argc, char **argv)
{
	if(argc < 4) {
		printf("Usage: ./FFEA_convert_traj_to_pdb [INPUT FFEA_TRAJECTORY FILE] [OUTPUT PDB FILE] [num_frames_to_convert] [FFEA scale (optional)]\n");
		return 1;
	}

	FILE *traj = NULL, *pdb_out = NULL;
	printf("Opening trajectory file %s for reading\n", argv[1]);
	if((traj = fopen(argv[1], "r")) == NULL) {
		printf("Error: Could not open trajectory file %s for reading.\n", argv[1]);
		return 1;
	}

	printf("Outputting to pdb file %s\n", argv[2]);
	if((pdb_out = fopen(argv[2], "w")) == NULL) {
		printf("Error: Could not open pdb file %s for reading.\n", argv[2]);
		return 1;
	}

	int num_frames_to_convert = atoi(argv[3]);
	double scale = 1.0;
	if(argc == 5) {
		scale /= atof(argv[4]);
	}

	int i, j;
	long long step = 0;
	int frame = 0, start = 0;;
	int **num_nodes, *num_conformations, num_blobs = 0;
	double x = 0, y = 0, z = 0;
	char line[255];
	char motion_state[7];

	// Scan through initial crap
	fscanf(traj, "%s\n", line);
	if(strcmp(line, "FFEA_trajectory_file") != 0) {
		printf("Error. Expected 'FFEA_trajectory_file' but got %s. May not be an FFEA_trajectory_file.\n", line);
		exit(0);
	}
	fscanf(traj, "\nInitialisation:\n");
	if(fscanf(traj, "Number of Blobs %d\n", &num_blobs) != 1) {
		printf("Error. Expected 'Number of Blobs x'\n");
		exit(0);
	}

	// Assign memory for conformations
	num_conformations = new int[num_blobs];
	num_nodes = new int*[num_blobs];
	char buf[100];
	fgets(buf, 25, traj);
	for(i = 0; i < num_blobs; ++i) {
		fgets(buf, 2, traj);
		num_conformations[i] = atoi(buf);
		num_nodes[i] = new int[num_conformations[i]];
	}
	fgets(buf,100, traj);
	for(i = 0; i < num_blobs; ++i) {
		if(fscanf(traj, "Blob %*d: Conformation %*d Nodes %d\n", &num_nodes[i][0]) != 1) {
			printf("Error. Expected 'Blob %%*d: Conformation %%*d Nodes %d\n", i);
			exit(0);
		}
		
	}
	fscanf(traj, "\n\n*\n");
	while(frame < num_frames_to_convert && !feof(traj)) {

		// Progress check
		frame++;
		if(frame % 100 == 0) {
			printf("step = %lld\n", step);
		}
		
		start = 0;
		fprintf(pdb_out, "MODEL     %4d\n", frame - 1);
		for(i = 0; i < num_blobs; ++i) {
			fgets(buf,255,traj);
			cout << buf << endl;
			exit(0);
			if(fscanf(traj, "Blob %*d, Conformation %*d, step %lld", &step) != 1) {
				printf("Error when reading 'Blob x, Conformation y, step z' line after or on step %lld\n", step);
				return -1;
			}

			// Rest of line
			fgets(buf, 255, traj);
			fscanf(traj, "%s\n", motion_state);
			if(strcmp(motion_state, "STATIC") == 0) {
				continue;
			}
			for(j = start; j < start + num_nodes[i][0]; j++) {

				// Read the x, y and z position values
				if(fscanf(traj, "%le %le %le", &x, &y, &z) != 3) {
					printf("Could not read coordinates on line %d of step %lld\n", (j+1), step);
					return -1;
				}

				// Skip the rest of the line
				while(fgetc(traj) != '\n');

				// FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)

				x *= scale;
				y *= scale;
				z *= scale;

				fprintf(pdb_out, "ATOM%7d  N   GLY     1    %8.3f%8.3f%8.3f\n", j, x, y, z);
			}
			start += num_nodes[i][0];
		}
		fscanf(traj, "*\nConformation Changes:\nBlob %*d: Conformation %*d -> Conformation %*d\n*\n");
		fprintf(pdb_out, "TER\n");
		fprintf(pdb_out, "ENDMDL\n");
	}

	fclose(traj);
	fclose(pdb_out);

	return 0;
}
