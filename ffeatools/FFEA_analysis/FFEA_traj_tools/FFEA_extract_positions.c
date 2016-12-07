// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	if(argc != 3) {
		printf("Usage: ./FFEA_extract_positions [INPUT TRAJECTORY FNAME] [XYZ OUTPUT FNAME]\n");
		return -1;
	}

	FILE *traj = NULL;
	FILE *out = NULL;
	int i;

	printf("traj = %s\n", argv[1]);
	printf("out = %s\n", argv[2]);

	if((traj = fopen(argv[1], "r")) == NULL) {
		printf("Could not open %s for reading\n", argv[1]);
	}

	if((out = fopen(argv[2], "w")) == NULL) {
		printf("Could not open %s for writing\n", argv[4]);
	}

	char nodestr[100];
	int num_nodes;
	int blob_id = -1;
	long long int step;
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
		i = 0;
		for(i = 0; i < num_nodes; i++) {
			if((result = fscanf(traj, "%lf %lf %lf %*f %*f %*f %*f %*f %*f %*f\n",&x, &y, &z)) != 3) {
					printf("result = %d\n", result);
					printf("Error reading node %d info at Blob %d, step %lld\n", i, blob_id, step);
					return -1;
			}
			fprintf(out, "%le %le %le\n", x, y, z);
		}
	}

	return 0;
}
