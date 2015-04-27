#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NEW_LINE_EVERY_10 if(num_on_line == 9) { fprintf(outx, "\n"); num_on_line = 0;	} else { num_on_line++;	}


int main(int argc, char **argv)
{
	if(argc < 4) {
		printf("Usage: ./FFEA_convert_traj_to_mdcrd [INPUT FFEA_TRAJECTORY FILE] [OUTPUT MDCRD FILE] [num frames to convert] [FFEA scale (optional)]\n");
		return 1;
	}

	FILE *traj = NULL, *outx = NULL;
	printf("Opening trajectory file %s for reading\n", argv[1]);
	if((traj = fopen(argv[1], "r")) == NULL) {
		printf("Error: Could not open trajectory file %s for reading.\n", argv[1]);
		return 1;
	}

	printf("Outputting to mdcrd file %s\n", argv[2]);
	if((outx = fopen(argv[2], "w")) == NULL) {
		printf("Error: Could not open mdcrd file %s for writing.\n", argv[2]);
		return 1;
	}

	fprintf(outx, "trajectory not created by ptraj\n");

	int num_frames_to_convert = atoi(argv[3]);
	double scale = 1.0;
	if(argc == 5) {
		scale /= atof(argv[4]);
	}
	int i, j;
	long long step = 0;
	int frame = 0;
	int *num_nodes, num_blobs = 0;
	double x = 0, y = 0, z = 0;
	int num_on_line = 0;
	char temp[10];
	char line[255];
	char motion_state[7];
//	int frame = 0;

	// Scan through initial crap
	fscanf(traj, "%s\n", &line);
	if(strcmp(line, "FFEA_trajectory_file") != 0) {
		printf("Error. Expected 'FFEA_trajectory_file' but got %s. May not be an FFEA_trajectory_file.\n", line);
		exit(0);
	}
	fscanf(traj, "\nInitialisation:\n");
	if(fscanf(traj, "Number of Blobs %d\n", &num_blobs) != 1) {
		printf("Error. Expected 'Number of Blobs x'\n");
		exit(0);
	}

	num_nodes = new int[num_blobs];

	for(i = 0; i < num_blobs; ++i) {
		if(fscanf(traj, "Blob %*d Nodes %d ", &num_nodes[i]) != 1) {
			printf("Error. Expected 'Blob %d Nodes ? \n", i);
			exit(0);
		}
		
	}
	fscanf(traj, "\n\n*\n");
	while(frame < num_frames_to_convert && !feof(traj)) {

		// Progress check
		frame++;
		if(frame % 1000 == 0) {
			printf("step = %lld\n", step);
		}

		for(i = 0; i < num_blobs; ++i) {
			if(fscanf(traj, "Blob %*d, Conformation %*d, step %lld\n", &step) != 1) {
				printf("Error when reading 'Blob x, Conformation y, step z' line after or on step %lld\n", step);
				return -1;
			}


			fscanf(traj, "%s\n", motion_state);
			if(strcmp(motion_state, "STATIC") == 0) {
				continue;
			}
			for(j = 0; j < num_nodes[i]; j++) {

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

				fprintf(outx, "%8.3f", x);
				NEW_LINE_EVERY_10
				fprintf(outx, "%8.3f", y);
				NEW_LINE_EVERY_10
				fprintf(outx, "%8.3f", z);
				NEW_LINE_EVERY_10
			}
			fprintf(outx, "\n");
			num_on_line = 0;
		}
		fscanf(traj, "*\n");
	}
	fclose(traj);
	fclose(outx);
	return 0;
}
