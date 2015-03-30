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

	if(argc != 6) {
		printf("Usage: make_cubic_mesh [num_nodes_x] [num_nodes_y] [num_nodes_z] [OUTPUT FILE NODES] [OUTPUT FILE TOPOLOGY]\n");
		return -1;
	}

	num_nodes_x = atoi(argv[1]);
	num_nodes_y = atoi(argv[2]);
	num_nodes_z = atoi(argv[3]);

	printf("Creating mesh with dimensions of %dx%dx%d nodes\n", num_nodes_x, num_nodes_y, num_nodes_z);

	FILE *nodes;
	if((nodes = fopen(argv[4], "w")) == NULL) {
		printf("Unable to open nodes file %s for writing:\n", argv[4]);
		perror(NULL);
		return -1;
	}

	fprintf(nodes, "walrus node file\n");
	fprintf(nodes, "num_nodes %d\n", num_nodes_x * num_nodes_y * num_nodes_z);
	fprintf(nodes, "num_surface_nodes ?\nnum_interior_nodes ?\nsurface nodes:\ninterior nodes:\n");

	// Create the nodes in a 3d grid
	for(k = 0; k < num_nodes_z; k++)
		for(j = 0; j < num_nodes_y; j++)
			for(i = 0; i < num_nodes_x; i++)
				fprintf(nodes, "%d %d %d\n", i, j, k);
	fclose(nodes);

	FILE *topology;
	if((topology = fopen(argv[5], "w")) == NULL) {
		printf("Unable to open topology file %s for writing:\n", argv[5]);
		perror(NULL);
		return -1;
	}


	fprintf(topology, "walrus topology file\n");
	fprintf(topology, "num_elements %d\n", (num_nodes_x - 1) * (num_nodes_y - 1) * (num_nodes_z - 1) * 5);
	fprintf(topology, "num_surface_elements ?\nnum_interior_elements ?\nsurface elements:\ninterior elements:\n");

	// Create the elements
	for(k = 0; k < num_nodes_z - 1; k++)
		for(j = 0; j < num_nodes_y - 1; j++)
			for(i = 0; i < num_nodes_x - 1; i++) {

				if(e_or_o(i)*e_or_o(j)*e_or_o(k) == 1) {
					even = 1;
					odd = 0;
				} else {
					even = 0;
					odd = 1;
				}

				// central tetrahedron
				fprintf(topology, "%d %d %d %d\n", 	NODE(i+odd,j+odd,k+odd),
									NODE(i+even,j+even,k+odd),
									NODE(i+even,j+odd,k+even),
									NODE(i+odd,j+even,k+even));

				// corner tetrahedrons
				fprintf(topology, "%d %d %d %d\n",	NODE(i+odd,j+odd,k+odd),
									NODE(i+even,j+odd,k+odd),
									NODE(i+even,j+even,k+odd),
									NODE(i+even,j+odd,k+even));

				fprintf(topology, "%d %d %d %d\n",	NODE(i+odd,j+odd,k+odd),
									NODE(i+odd,j+odd,k+even),
									NODE(i+odd,j+even,k+even),
									NODE(i+even,j+odd,k+even));

				fprintf(topology, "%d %d %d %d\n",	NODE(i+even,j+even,k+odd),
									NODE(i+even,j+odd,k+even),
									NODE(i+odd,j+even,k+even),
									NODE(i+even,j+even,k+even));

				fprintf(topology, "%d %d %d %d\n",	NODE(i+odd,j+odd,k+odd),
									NODE(i+odd,j+even,k+even),
									NODE(i+even,j+even,k+odd),
									NODE(i+odd,j+even,k+odd));
		}

	fclose(topology);
}
