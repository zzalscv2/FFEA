#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int count_commas(char *s)
{
	int i = 0;
	int count = 0;
	for(;;) {
		if(s[i] == ',') {
			count += 1;
		} else if(s[i] == '\0') {
			break;
		}
		i += 1;
	}

	return count;
}

int main(int argc, char **argv)
{
	if(argc != 5) {
		printf("Usage: ./extract [INPUT TRAJECTORY FNAME] [BLOB NUMBER] [NODE NUMBERS (COMMA SEPARATED)] [OUTPUT TRAJECTORY FNAME]\n");
		return -1;
	}

	FILE *traj = NULL;
	FILE *out = NULL;
	int blob_num = -1;
	int i;

	blob_num = atoi(argv[2]);

	int num_extract_nodes = count_commas(argv[3]) + 1;

	printf("num_extract_nodes = %d\n", num_extract_nodes);

	int *nodes_list = malloc(num_extract_nodes * sizeof(int));

	int ni = 0;
	int isi = 0;
	char int_str[256];
	i = 0;
	for(;;) {
		if (argv[3][i] == '\0') {
			int_str[isi] = '\0';
			nodes_list[ni] = atoi(int_str);
			break;

		} else if(argv[3][i] == ',') {
			int_str[isi] = '\0';
			nodes_list[ni] = atoi(int_str);
			ni += 1;
			isi = 0;

			if(ni == num_extract_nodes) {
				break;
			}
		} else {
			int_str[isi] = argv[3][i];
			isi += 1;
		}
		i += 1;
	}

	printf("Nodes to extract:\n");
	for(i = 0; i < num_extract_nodes; i++) {
		printf("%d\n", nodes_list[i]);
		if(i != 0) {
			if(nodes_list[i] < nodes_list[i-1]) {
				printf("Error: Node indices must be in ascending order.\n");
				return 0;
			} else if(nodes_list[i] == nodes_list[i-1]) {
				printf("Error: Duplicate node indicies.\n");
				return 0;
			}
		}
	}

	printf("traj = %s\n", argv[1]);
	printf("blob_num = %d\n", blob_num);
	printf("out = %s\n", argv[4]);

	if((traj = fopen(argv[1], "r")) == NULL) {
		printf("Could not open %s for reading\n", argv[1]);
	}

	if((out = fopen(argv[4], "w")) == NULL) {
		printf("Could not open %s for writing\n", argv[4]);
	}

	int blob_id = -1;
	long long int step = -1;
	char nodestr[100];
	int num_nodes;
	double x = 0, y = 0, z = 0;
	char c;
	int result;
	int current;
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
		fprintf(out, "*\nBlob %d, step %lld\n%d\n", blob_id, step, num_extract_nodes);
		if(blob_id == blob_num) {
			i = 0;
			current = 0;
			while(i < num_nodes) {
				if((result = fscanf(traj, "%lf %lf %lf %*f %*f %*f %*f %*f %*f %*f\n",&x, &y, &z)) != 3) {
					printf("result = %d\n", result);
					printf("Error reading node %d info at Blob %d, step %lld\n", i, blob_id, step);
					return -1;
				}
				if(current < num_extract_nodes) {
					if(i == nodes_list[current]) {
						fprintf(out, "%le %le %le 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00\n", x, y, z);
						current += 1;
					}
				}
				i += 1;
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
