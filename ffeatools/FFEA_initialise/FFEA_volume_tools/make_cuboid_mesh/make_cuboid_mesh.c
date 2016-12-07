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

#include <stdlib.h>
#include <stdio.h>

#define NODE(I,J,K) ((I) + (J) * num_nodes_x + (K) * num_nodes_x * num_nodes_y)

static inline int e_or_o(int i)
{
	if(i%2 == 0) return 1;
	else return -1;
}

int main(int argc, char *argv[])
{
	int i, j, k, even, odd;
	int num_nodes_x = 0, num_nodes_y = 0, num_nodes_z = 0;

	if(argc != 5) {
		printf("Usage: make_cuboid_mesh [num_nodes_x] [num_nodes_y] [num_nodes_z] [OUTPUT .DAT MESH FNAME]\n");
		return -1;
	}

	num_nodes_x = atoi(argv[1]);
	num_nodes_y = atoi(argv[2]);
	num_nodes_z = atoi(argv[3]);

	printf("Creating .dat mesh with dimensions of %dx%dx%d nodes\n", num_nodes_x, num_nodes_y, num_nodes_z);

	FILE *datfile;
	if((datfile = fopen(argv[4], "w")) == NULL) {
		printf("Unable to open datfile file %s for writing:\n", argv[4]);
		perror(NULL);
		return -1;
	}

	int num_nodes = num_nodes_x * num_nodes_y * num_nodes_z;
	int num_elements = (num_nodes_x - 1) * (num_nodes_y - 1) * (num_nodes_z - 1) * 5;
	fprintf(datfile, "%d %d\n", num_nodes, num_elements);

	// Create the nodes in a 3d grid
	int node_index = 1;
	for(k = 0; k < num_nodes_z; k++)
		for(j = 0; j < num_nodes_y; j++)
			for(i = 0; i < num_nodes_x; i++)
				fprintf(datfile, "%d %d %d %d\n", node_index, i, j, k);
				node_index += 1;

	// Create the elements
	int element_index = 1;
	for(k = 0; k < num_nodes_z - 1; k++) {
		for(j = 0; j < num_nodes_y - 1; j++) {
			for(i = 0; i < num_nodes_x - 1; i++) {

				if(e_or_o(i)*e_or_o(j)*e_or_o(k) == 1) {
					even = 1;
					odd = 0;
				} else {
					even = 0;
					odd = 1;
				}

				// central tetrahedron
				fprintf(datfile, "%d %d %d %d %d\n",	element_index,
									NODE(i+odd,j+odd,k+odd) + 1,
									NODE(i+even,j+even,k+odd) + 1,
									NODE(i+even,j+odd,k+even) + 1,
									NODE(i+odd,j+even,k+even) + 1);
				element_index += 1;

				// corner tetrahedrons
				fprintf(datfile, "%d %d %d %d %d\n",	element_index,
									NODE(i+odd,j+odd,k+odd) + 1,
									NODE(i+even,j+odd,k+odd) + 1,
									NODE(i+even,j+even,k+odd) + 1,
									NODE(i+even,j+odd,k+even) + 1);
				element_index += 1;


				fprintf(datfile, "%d %d %d %d %d\n",	element_index,
									NODE(i+odd,j+odd,k+odd) + 1,
									NODE(i+odd,j+odd,k+even) + 1,
									NODE(i+odd,j+even,k+even) + 1,
									NODE(i+even,j+odd,k+even) + 1);
				element_index += 1;

				fprintf(datfile, "%d %d %d %d %d\n",	element_index,
									NODE(i+even,j+even,k+odd) + 1,
									NODE(i+even,j+odd,k+even) + 1,
									NODE(i+odd,j+even,k+even) + 1,
									NODE(i+even,j+even,k+even) + 1);
				element_index += 1;

				fprintf(datfile, "%d %d %d %d %d\n",	element_index,
									NODE(i+odd,j+odd,k+odd) + 1,
									NODE(i+odd,j+even,k+even) + 1,
									NODE(i+even,j+even,k+odd) + 1,
									NODE(i+odd,j+even,k+odd) + 1);
				element_index += 1;
			}
		}
	}

	fclose(datfile);
}
