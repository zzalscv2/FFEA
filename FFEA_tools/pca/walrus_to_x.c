#include <stdio.h>
#include <stdlib.h>

#define SCALE 1e9
#define NEW_LINE_EVERY_10 if(num_on_line == 9) { fprintf(outx, "\n"); num_on_line = 0;	} else { num_on_line++;	}


int main(int argc, char **argv)
{
	if(argc != 3) {
		printf("Usage: ./walrus_to_x [INPUT WALRUS_TRAJECTORY FILE] [OUTPUT PDB FILE]\n");
		return 1;
	}

	FILE *traj = NULL, *outx = NULL;
	printf("Opening trajectory file %s for reading\n", argv[1]);
	if((traj = fopen(argv[1], "r")) == NULL) {
		printf("Error: Could not open trajectory file %s for reading.\n", argv[1]);
		return 1;
	}

	printf("Outputting to x file %s\n", argv[2]);
	if((outx = fopen(argv[2], "w")) == NULL) {
		printf("Error: Could not open trajectory file %s for reading.\n", argv[2]);
		return 1;
	}

	fprintf(outx, "trajectory not created by ptraj\n");

	long long step = 0;
	int num_nodes = 0;
	double x = 0, y = 0, z = 0;
	int num_on_line = 0;
	int every_1000 = 0;
	char temp[10];
//	int frame = 0;
	while(!feof(traj)) {
		if(fscanf(traj, "*\nBlob %*d, step %lld\n", &step) != 1) {
			printf("Error when reading '*\\nBlob x, step y\\n' line after or on step %lld\n", step);
			break;
			return -1;
		}
		every_1000++;
		if(every_1000 == 1000) {
			printf("step = %lld\n", step);
			every_1000 = 0;
		}
		// STATIC etc
		fscanf(traj, "%s\n", &temp);
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

			fprintf(outx, "%8.3f", x);
			NEW_LINE_EVERY_10
			fprintf(outx, "%8.3f", y);
			NEW_LINE_EVERY_10
			fprintf(outx, "%8.3f", z);
			NEW_LINE_EVERY_10
		}
		fprintf(outx, "\n");
		num_on_line = 0;

//		frame++;
//		if(frame > 300) {
//			break;
//		}
	}

	fclose(traj);
	fclose(outx);

	return 0;
}
