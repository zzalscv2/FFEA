#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	if(argc != 5) {
		printf("Usage: ./extract [INPUT TRAJECTORY FNAME] [BLOB NUMBER] [NODE NUMBER] [OUTPUT POSITIONS FNAME]\n");
		return -1;
	}

	FILE *traj = NULL;
	FILE *out = NULL;
	int blob_num = -1;
	int node_num = -1;

	blob_num = atoi(argv[2]);
	node_num = atoi(argv[3]);

	printf("traj = %s\n", argv[1]);
	printf("blob_num = %d\n", blob_num);
	printf("node_num = %d\n", node_num);
	printf("out = %s\n", argv[4]);

	if((traj = fopen(argv[1], "r")) == NULL) {
		printf("Could not open %s for reading\n", argv[1]);
	}

	if((out = fopen(argv[4], "w")) == NULL) {
		printf("Could not open %s for writing\n", argv[4]);
	}

	int i;
	int blob_id = -1;
	long long int step = -1;
	char nodestr[100];
	int num_nodes;
	double x = 0, y = 0, z = 0;
	char c;
	int result;
	for(;;) {
		if((c = fgetc(traj)) != '*') {
			ungetc(c, traj);
		} else {
			fgetc(traj);
		}
		if(feof(traj)) {
			printf("Reached end of file. No errors.");
			break;
		}
		if(fscanf(traj, "Blob %d, step %lld\n", &blob_id, &step) != 2) {
			printf("Error reading header info after Blob %d at step %lld\n", blob_id, step);
			return -1;
		}
		if(fgets(nodestr, 100, traj) == NULL) {
			printf("Error\n");
			return -1;
		}
		if(strcmp(nodestr, "STATIC\n") == 0) {
			continue;
		}
		num_nodes = atoi(nodestr);
		if(blob_id == blob_num) {
			for(i = 0; i < node_num; i++) {
				if((result = fscanf(traj, "%lf %lf %lf %*f %*f %*f %*f %*f %*f %*f\n",&x, &y, &z)) != 3) {
					printf("result = %d\n", result);
					printf("Error reading node %d info at Blob %d, step %lld\n", i, blob_id, step);
					return -1;
				}
			}
			if((result = fscanf(traj, "%lf %lf %lf %*f %*f %*f %*f %*f %*f %*f\n",&x, &y, &z)) != 3) {
				printf("result = %d\n", result);
				printf("Error reading node %d info at Blob %d, step %lld\n", i, blob_id, step);
				return -1;
			}
			fprintf(out, "%le %le %le\n", x, y, z);
			for(i = node_num + 1; i < num_nodes; i++) {
				if((result = fscanf(traj, "%lf %lf %lf %*f %*f %*f %*f %*f %*f %*f\n",&x, &y, &z)) != 3) {
					printf("result = %d\n", result);
					printf("Error reading node %d info at Blob %d, step %lld\n", i, blob_id, step);
					return -1;
				}
			}

		} else {
			for(i = 0; i < num_nodes; i++) {
				if((result = fscanf(traj, "%lf %lf %lf %*f %*f %*f %*f %*f %*f %*f\n",&x, &y, &z)) != 3) {
					printf("result = %d\n", result);
					printf("Error reading node %d info at Blob %d, step %lld\n", i, blob_id, step);
					return -1;
				}
			}
		}
	}

	return 0;
}
