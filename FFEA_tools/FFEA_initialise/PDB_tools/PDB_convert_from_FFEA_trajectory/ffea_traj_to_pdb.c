#include <stdio.h>
#include <stdlib.h>

#define SCALE 1e10

int main(int argc, char **argv)
{
	if(argc != 3) {
		printf("Usage: ./ffea_traj_to_pdb [INPUT WALRUS_TRAJECTORY FILE] [OUTPUT PDB FILE]\n");
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

	long long step = 0;
	int num_nodes = 0;
	double x = 0, y = 0, z = 0;
	while(!feof(traj)) {
		if(fscanf(traj, "*\nBlob %*d, step %lld\n", &step) != 1) {
			printf("Error when reading '*\\nBlob x, step y\\n' line after or on step %lld\n", step);
			return -1;
		}

		if(fscanf(traj, "%d\n", &num_nodes) != 1) {
			printf("Could not read number of nodes line at step %lld\n", step);
			return -1;
		}

		for(int i = 0; i < num_nodes; i++) {
			// Read the x, y and z position values
			if(fscanf(traj, "%le %le %le", &x, &y, &z) != 3) {
				printf("Could not read coordinates on line %d of step %lld\n", (i+1), step);
				return -1;
			}

			// Skip the rest of the line
			while(fgetc(traj) != '\n');

			// FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)

			x *= SCALE;
			y *= SCALE;
			z *= SCALE;

			fprintf(pdb_out, "ATOM%7d  N   GLY     1    %8.3f%8.3f%8.3f\n", i, x, y, z);
		}
		fprintf(pdb_out, "TER\nENDMDL\n");
	}

	fclose(traj);
	fclose(pdb_out);

	return 0;
}
