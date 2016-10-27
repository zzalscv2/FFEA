#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "marching_cube_lookup.h"
#include "euler_characteristic.h"

#define MAX_FNAME_LENGTH 255

#define MAP_INDEX(X, Y, Z) (((Z)*map->NY*map->NX) + ((Y)*map->NX) + (X))

// As defined by the spec
#define HEADER_SIZE_IN_WORDS 256
#define WORD_SIZE 4

#define OUTPUT_FORMAT_TEXT	0
#define OUTPUT_FORMAT_SURF	1
#define OUTPUT_FORMAT_OBJ	2
#define OUTPUT_FORMAT_OFF	3
#define OUTPUT_FORMAT_STL	4
#define OUTPUT_FORMAT_MAP	5

// Error text in colour
void Error_text() {
    printf("\e[31mERROR: \e[m");
}

// Hold the header info
typedef struct
{
	signed int NC, NR, NS;
	signed int MODE;
	signed int NCSTART, NRSTART, NSSTART;
	signed int NX, NY, NZ;
	float X_LENGTH, Y_LENGTH, Z_LENGTH;
	float ALPHA, BETA, GAMMA;
	signed int MAPC, MAPR, MAPS;
	float AMIN, AMAX, AMEAN;
	signed int ISPG;
	signed int NSYMBT;
	signed int LSKFLG;
	float SKWMAT[9];
	float SKWTRN[3];
	uint32_t EXTRA[15];
	char MAP[5];
	uint32_t MACHST;
	float RMS;
	signed int NLABL;
	char LABEL_N[200 * WORD_SIZE + 1];
} header_info;

// Hold the map data
typedef struct
{
	int num_voxels;
	int NX, NY, NZ;
	float voxel_sizeX, voxel_sizeY, voxel_sizeZ;
	float *data;
} map_data;

class vector3
{
	public:
		double x, y, z;

		void scale(float s, float t, float u)
		{
			x *= s;
			y *= t;
			z *= u;
		}

		vector3 subtract(vector3 a) {

			vector3 s;
			s.x = x - a.x;
			s.y = y - a.y;
			s.z = z - a.z;

			return s;	
		}

		vector3 cross(vector3 a) {
		
			vector3 c;
			c.x = y * a.z - z * a.y;
			c.y = z * a.x - x * a.z;
			c.z = x * a.y - y * a.x;

			return c;
		}

		void normalise() {
			double mag = 0.0;
			
			mag = sqrt(x * x + y * y + z * z);
			x /= mag;
			y /= mag;
			z /= mag;	
		}

};

class surf_face
{
	public:
		surf_face(int n1, int n2, int n3)
		{
			this->n1 = n1;
			this->n2 = n2;
			this->n3 = n3;
		}

		vector3 get_normal(std::vector<vector3*> nodes) {
			
			vector3 n, e1, e2;

			// Get edge vectors
			e1 = nodes[n2]->subtract(*nodes[n1]);
			e2 = nodes[n3]->subtract(*nodes[n1]);

			// Cross for normal
			n = e1.cross(e2);

			n.normalise();
			return n;
		}
		int n1, n2, n3;
};


class grid_edge
{
	public:
		grid_edge()
		{
			voxel_index_1 = 0;
			voxel_index_2 = 0;
		}

		void set_edge(int index1, int index2)
		{
			if(index1 < index2) {
				voxel_index_1 = index1;
				voxel_index_2 = index2;
			} else {
				voxel_index_1 = index2;
				voxel_index_2 = index1;
			}
		}

		void set_edge(grid_edge *e)
		{
			voxel_index_1 = e->get_voxel_index_1();
			voxel_index_2 = e->get_voxel_index_2();
		}

		int get_voxel_index_1()
		{
			return voxel_index_1;
		}

		int get_voxel_index_2()
		{
			return voxel_index_2;
		}

		bool equals(grid_edge *e)
		{
			if(voxel_index_1 != e->get_voxel_index_1()) {
				return false;
			}

			if(voxel_index_2 != e->get_voxel_index_2()) {
				return false;
			}

			return true;
		}

	private:
		int voxel_index_1, voxel_index_2;
};

class edge_node_lookup
{
	public:
		edge_node_lookup(grid_edge *e, int node_index)
		{
			edge.set_edge(e);
			this->node_index = node_index;
		}

		bool equals(grid_edge *e)
		{
			return edge.equals(e);
		}

		bool equals_voxel_index_2(int vi2)
		{
			if(edge.get_voxel_index_2() == vi2) {
				return true;
			} else {
				return false;
			}
		}

		int get_node_index()
		{
			return node_index;
		}

	private:
		grid_edge edge;
		int node_index;
};

class map_to_surf_generator
{
	public:
		map_to_surf_generator() {
			num_nodes = 0;
			num_faces = 0;
		}

		void generate_surface(map_data *map, float level, int interpolate) {
			this->map = map;
			edge_node_lookup_table = new std::vector<edge_node_lookup*>[map->num_voxels];

			//           __________
			//          /|        /|
			//         /_________/ |
			//         | |       | |
			//         | |_______|_|
			//         | /       | /
			//         |/________|/
			//

			int voxel_index[8];

			for(int z = 0; z < map->NZ - 1; z++) {
				for(int y = 0; y < map->NY - 1; y++) {
					for(int x = 0; x < map->NX - 1; x++) {
						voxel_index[0] = MAP_INDEX(x    , y,     z    );
						voxel_index[1] = MAP_INDEX(x + 1, y,     z    );
						voxel_index[2] = MAP_INDEX(x + 1, y + 1, z    );
						voxel_index[3] = MAP_INDEX(x    , y + 1, z    );
						voxel_index[4] = MAP_INDEX(x,     y,     z + 1);
						voxel_index[5] = MAP_INDEX(x + 1, y,     z + 1);
						voxel_index[6] = MAP_INDEX(x + 1, y + 1, z + 1);
						voxel_index[7] = MAP_INDEX(x    , y + 1, z + 1);

						generate_faces_for_cube(voxel_index, level, interpolate);
					}
				}
				printf("Generating surface: %.0f%%\n", (z/(float)(map->NZ - 1))*100);
				printf("\033[F\033[J"); // Replace previous "progress" line
			}
			printf("Generating surface: 100%%\n");

			printf("Scaling by voxel size (%f x %f x %f Angstroms)\n", map->voxel_sizeX, map->voxel_sizeY, map->voxel_sizeZ);
			for(int i = 0; i < num_nodes; i++) {
				nodes[i]->scale(map->voxel_sizeX, map->voxel_sizeY, map->voxel_sizeZ);
			}
		}

		void write_netgen_surf_file(char *output_surface_fname, double cent_x, double cent_y, double cent_z)
		{
			FILE *surffile = NULL;
			if((surffile = fopen(output_surface_fname, "w")) == NULL) {
				printf("Could not open output file %s for writing\n", output_surface_fname);
				return;
			}

			fprintf(surffile, "surfacemesh\n");
			fprintf(surffile, "%d\n", num_nodes);
			for(int i = 0; i < num_nodes; i++) {
				fprintf(surffile, "%e %e %e\n", nodes[i]->x + cent_x, nodes[i]->y + cent_y, nodes[i]->z + cent_z);
			}

			fprintf(surffile, "%d\n", num_faces);
			for(int i = 0; i < num_faces; i++) {
				fprintf(surffile, "%d %d %d\n", faces[i]->n1 + 1, faces[i]->n2 + 1, faces[i]->n3 + 1);
			}

			fclose(surffile);

			printf("Wrote %d nodes and %d faces in netgen surf format to %s\n", num_nodes, num_faces, output_surface_fname);
		}

		void write_OFF_file(char *output_surface_fname, double cent_x, double cent_y, double cent_z)
		{
			FILE *surffile = NULL;
			if((surffile = fopen(output_surface_fname, "w")) == NULL) {
				printf("Could not open output file %s for writing\n", output_surface_fname);
				return;
			}

			fprintf(surffile, "OFF\n");
			fprintf(surffile, "%d %d 0\n", num_nodes, num_faces);
			for(int i = 0; i < num_nodes; i++) {
				fprintf(surffile, "%e %e %e\n", nodes[i]->x + cent_x, nodes[i]->y + cent_y, nodes[i]->z + cent_z);
			}
			for(int i = 0; i < num_faces; i++) {
				fprintf(surffile, "3 %d %d %d\n", faces[i]->n1, faces[i]->n2, faces[i]->n3);
			}
			fclose(surffile);

			printf("Wrote %d nodes and %d faces in OFF format to %s\n", num_nodes, num_faces, output_surface_fname);
		}

		void write_STL_file(char *output_surface_fname, double cent_x, double cent_y, double cent_z)
		{
			FILE *surffile = NULL;
			if((surffile = fopen(output_surface_fname, "w")) == NULL) {
				printf("Could not open output file %s for writing\n", output_surface_fname);
				return;
			}

			fprintf(surffile, "solid ffeatools_STL\n\n");
			vector3 n;
			for(int i = 0; i < num_faces; i++) {

				// Get normal first
				n = faces[i]->get_normal(nodes);
				fprintf(surffile, "facet normal %3.2f %3.2f %3.2f\n", n.x, n.y, n.z);
				fprintf(surffile, "\touter loop\n");
				fprintf(surffile, "\t\tvertex %3.6f %3.6f %3.6f\n", nodes[faces[i]->n1]->x + cent_x, nodes[faces[i]->n1]->y + cent_y, nodes[faces[i]->n1]->z + cent_z);
				fprintf(surffile, "\t\tvertex %3.6f %3.6f %3.6f\n", nodes[faces[i]->n2]->x + cent_x, nodes[faces[i]->n2]->y + cent_y, nodes[faces[i]->n2]->z + cent_z);
				fprintf(surffile, "\t\tvertex %3.6f %3.6f %3.6f\n", nodes[faces[i]->n3]->x + cent_x, nodes[faces[i]->n3]->y + cent_y, nodes[faces[i]->n3]->z + cent_z);
				fprintf(surffile, "\tendloop\n");
				fprintf(surffile, "endfacet\n\n");
			}

			fprintf(surffile, "endsolid ffeatools_STL\n");
				//fprintf(surffile, "f %d//%d %d//%d %d//%d\n", faces[i]->n1 + 1, faces[i]->n1 + 1, faces[i]->n2 + 1, faces[i]->n2 + 1, faces[i]->n3 + 1, faces[i]->n3 + 1);
			fclose(surffile);

			printf("Wrote %d nodes and %d faces in STL format to %s\n", num_nodes, num_faces, output_surface_fname);
		}

		void write_OBJ_file(char *output_surface_fname, double cent_x, double cent_y, double cent_z)
		{
			FILE *surffile = NULL;
			if((surffile = fopen(output_surface_fname, "w")) == NULL) {
				printf("Could not open output file %s for writing\n", output_surface_fname);
				return;
			}

			fprintf(surffile, "# Created by emdb_map_to_ffea.c\n");
			fprintf(surffile, "# Number of vertices: %d\n", num_nodes);
			for(int i = 0; i < num_nodes; i++) {
				fprintf(surffile, "v %e %e %e\n", nodes[i]->x + cent_x, nodes[i]->y + cent_y, nodes[i]->z + cent_z);
			}
			fprintf(surffile, "# Number of triangles: %d\n", num_faces);
			for(int i = 0; i < num_faces; i++) {
				fprintf(surffile, "f %d//%d %d//%d %d//%d\n", faces[i]->n1 + 1, faces[i]->n1 + 1, faces[i]->n2 + 1, faces[i]->n2 + 1, faces[i]->n3 + 1, faces[i]->n3 + 1);
			}

			fclose(surffile);

			printf("Wrote %d nodes and %d faces in OBJ format to %s\n", num_nodes, num_faces, output_surface_fname);
		}

	private:
		map_data *map;

		int num_nodes, num_faces;

		std::vector<vector3*> nodes;
		std::vector<surf_face*> faces;

		std::vector<edge_node_lookup*> *edge_node_lookup_table;

		void add_tri(grid_edge *e1, grid_edge *e2, grid_edge *e3, float level, int interpolate)
		{
			int n1 = add_node_to_edge(e1, level, interpolate);
			int n2 = add_node_to_edge(e2, level, interpolate);
			int n3 = add_node_to_edge(e3, level, interpolate);
			surf_face *new_face = new surf_face(n1, n2, n3);
			faces.push_back(new_face);
			num_faces++;
		}

		void map_index_to_position(int index, vector3 *v)
		{
			int x, y, z;
			x = index % map->NX;
			y = ((index - x)/map->NX) % map->NY;
			z = (index - x - y * map->NX)/(map->NX * map->NY);

			v->x = x;
			v->y = y;
			v->z = z;
		}

		int add_node_to_edge(grid_edge *e, float level, int interpolate)
		{
			// search for existing node on this edge
//			for(int i = 0; i < num_nodes; i++) {
//				if(edge_nodes[i]->equals(e) == true) {
//					// if one exists, reuse it
//					return i;
//				}
//			}

			int vindex_1 = e->get_voxel_index_1(), vindex_2 = e->get_voxel_index_2();

			for(unsigned int i = 0; i < edge_node_lookup_table[vindex_1].size(); i++) {
				if((edge_node_lookup_table[vindex_1])[i]->equals_voxel_index_2(vindex_2) == true) {
					// if one exists, reuse it
					return (edge_node_lookup_table[vindex_1])[i]->get_node_index();
				}
			}

			// if no such node, create it:
			int new_node_index = num_nodes;

			// get the end voxel positions
			vector3 v1, v2;
			map_index_to_position(vindex_1, &v1);
			map_index_to_position(vindex_2, &v2);

			// get the new node position by linearly interpolating between these points
			vector3 *new_node_pos = new vector3();
			if(interpolate == 1) {
				calc_node_position(&v1, &v2, map->data[vindex_1], map->data[vindex_2], level, new_node_pos);
			} else {
				calc_node_position(&v1, &v2, new_node_pos);
			}

			// add new node to the node array
			nodes.push_back(new_node_pos);

			// create the edge-node lookup entry and add it to the lookup table
			edge_node_lookup *new_edge_node_lookup = new edge_node_lookup(e, new_node_index);
			edge_node_lookup_table[vindex_1].push_back(new_edge_node_lookup);

			num_nodes++;

			return new_node_index;
		}

		// Place the node at the midpoint between p1 and p2
		void calc_node_position(vector3 *p1, vector3 *p2, vector3 *return_v)
		{
			return_v->x = (p1->x + p2->x)/2;
			return_v->y = (p1->y + p2->y)/2;
			return_v->z = (p1->z + p2->z)/2;
		}

		// Linearly interpolate between p1 and p2 based on the density value at each point
		void calc_node_position(vector3 *p1, vector3 *p2, float density_1, float density_2, float level, vector3 *return_v)
		{
			// Avoid nasty NaNs...
			if(density_2 == -INFINITY) {
				return_v->x = p1->x;
				return_v->y = p1->y;
				return_v->z = p1->z;
				return;
			}
			if(density_1 == -INFINITY) {
				return_v->x = p2->x;
				return_v->y = p2->y;
				return_v->z = p2->z;
				return;
			}

			float s = (level - density_1) / (density_2 - density_1);

			return_v->x = p1->x + (p2->x - p1->x) * s;
			return_v->y = p1->y + (p2->y - p1->y) * s;
			return_v->z = p1->z + (p2->z - p1->z) * s;
		}

		// Generate the appropriate faces for the given cube of 8 voxels
		void generate_faces_for_cube(int voxel_index[8], float level, int interpolate)
		{
			int cube_case;
			grid_edge edge_id[12];

			// Get the "case" for this cube (A byte between 0 and 255 denoting which of the 2^8 different configurations this cube is in)
			cube_case = 0;
			for(int i = 0; i < 8; i++) {
				if(map->data[voxel_index[i]] <= level) cube_case |= (1 << i);
			}

			// If all voxels are filled or unfilled, do nothing
			if(edgeTable[cube_case] == 0) {
				return;
			}

			// Otherwise, look up the nodes at which the surface intersects this cube
			if(edgeTable[cube_case] & 1) {
				edge_id[0].set_edge(voxel_index[0], voxel_index[1]);
			}
			if(edgeTable[cube_case] & 2) {
				edge_id[1].set_edge(voxel_index[1], voxel_index[2]);
			}
			if(edgeTable[cube_case] & 4) {
				edge_id[2].set_edge(voxel_index[2], voxel_index[3]);
			}
			if(edgeTable[cube_case] & 8) {
				edge_id[3].set_edge(voxel_index[3], voxel_index[0]);
			}
			if(edgeTable[cube_case] & 16) {
				edge_id[4].set_edge(voxel_index[4], voxel_index[5]);
			}
			if(edgeTable[cube_case] & 32) {
				edge_id[5].set_edge(voxel_index[5], voxel_index[6]);
			}
			if(edgeTable[cube_case] & 64) {
				edge_id[6].set_edge(voxel_index[6], voxel_index[7]);
			}
			if(edgeTable[cube_case] & 128) {
				edge_id[7].set_edge(voxel_index[7], voxel_index[4]);
			}
			if(edgeTable[cube_case] & 256) {
				edge_id[8].set_edge(voxel_index[0], voxel_index[4]);
			}
			if(edgeTable[cube_case] & 512) {
				edge_id[9].set_edge(voxel_index[1], voxel_index[5]);
			}
			if(edgeTable[cube_case] & 1024) {
				edge_id[10].set_edge(voxel_index[2], voxel_index[6]);
			}
			if(edgeTable[cube_case] & 2048) {
				edge_id[11].set_edge(voxel_index[3], voxel_index[7]);
			}

			// Add the triangular faces
			int i = 0;
			while(triTable[cube_case][i] != -1) {
				add_tri(&edge_id[triTable[cube_case][i]], &edge_id[triTable[cube_case][i + 1]], &edge_id[triTable[cube_case][i + 2]], level, interpolate);
				i += 3;
			}
		}
};

class map_flood_filler
{
	public:
		int cull_floaters(map_data *map, float level, int cull_threshold) {
			int *workspace = new int[map->num_voxels];
			for(int i = 0; i < map->NX; i++) {
				for(int j = 0; j < map->NY; j++) {
					for(int k = 0; k < map->NZ; k++) {
						if(i > 1 && i < map->NX - 1 && j > 1 && j < map->NY - 1 && k > 1 && k < map->NZ - 1) {
							workspace[MAP_INDEX(i, j, k)] = 0;
						} else {
							// map border
							workspace[MAP_INDEX(i, j, k)] = -1;
						}
					}
				}
			}

			int num_culled = 0, num_kept = 0;
			for(int i = 1; i < map->NX - 1; i++) {
				for(int j = 1; j < map->NY - 1; j++) {
					for(int k = 1; k < map->NZ - 1; k++) {
						if(workspace[MAP_INDEX(i, j, k)] == 0) {
							if(map->data[MAP_INDEX(i,j,k)] > level) {
								if(flood_fill_and_cull(i, j, k, map, level, workspace, cull_threshold) == true) {
									num_culled++;
								} else {
									num_kept++;
								}
							}
						}
					}
				}
			}
			printf("Culled %d 'floaters' (size < %d voxels)\nKept %d.\n", num_culled, cull_threshold, num_kept);
			delete[] workspace;
			return 0;
		}

		int fill_cavities(map_data *map, float level, int fill_threshold) {
			int *workspace = new int[map->num_voxels];
			for(int i = 0; i < map->NX; i++) {
				for(int j = 0; j < map->NY; j++) {
					for(int k = 0; k < map->NZ; k++) {
						if(i > 1 && i < map->NX - 1 && j > 1 && j < map->NY - 1 && k > 1 && k < map->NZ - 1) {
							workspace[MAP_INDEX(i, j, k)] = 0;
						} else {
							// map border
							workspace[MAP_INDEX(i, j, k)] = -1;
						}
					}
				}
			}

			int num_filled = 0, num_left_empty = 0;
			for(int i = 1; i < map->NX - 1; i++) {
				for(int j = 1; j < map->NY - 1; j++) {
					for(int k = 1; k < map->NZ - 1; k++) {
						if(workspace[MAP_INDEX(i, j, k)] == 0) {
							if(map->data[MAP_INDEX(i,j,k)] <= level) {
								if(flood_fill_cavity(i, j, k, map, level, workspace, fill_threshold) == true) {
									num_filled++;
								} else {
									num_left_empty++;
								}
							}
						}
					}
				}
			}
			printf("Filled %d holes (size < %d voxels)\nKept %d.\n", num_filled, fill_threshold, num_left_empty);
			delete[] workspace;
			return 0;
		}

	private:

		class map_pos
		{
			public:
				map_pos(int i, int j, int k) {
					this->i = i;
					this->j = j;
					this->k = k;
				}
				int i, j, k;
		};

		bool flood_fill_and_cull(int i, int j, int k, map_data *map, float level, int *workspace, int cull_threshold) {
			std::vector<map_pos> to_check;
			std::vector<map_pos> to_delete;
			to_check.push_back(map_pos(i,j,k));
			int num_filled = 0;
			while(to_check.size() > 0) {
				map_pos p = to_check.back();
				if(workspace[MAP_INDEX(p.i, p.j, p.k)] != 0) {
					to_check.pop_back();
				} else {
					if(map->data[MAP_INDEX(p.i, p.j, p.k)] > level) {
						workspace[MAP_INDEX(p.i, p.j, p.k)] = 1;
						to_delete.push_back(map_pos(p.i, p.j, p.k));
						num_filled++;
						to_check.push_back(map_pos(p.i-1,p.j,p.k));
						to_check.push_back(map_pos(p.i+1,p.j,p.k));
						to_check.push_back(map_pos(p.i,p.j-1,p.k));
						to_check.push_back(map_pos(p.i,p.j+1,p.k));
						to_check.push_back(map_pos(p.i,p.j,p.k-1));
						to_check.push_back(map_pos(p.i,p.j,p.k+1));
					} else {
						to_check.pop_back();
					}
				}
			}

			if(num_filled <= cull_threshold) {
				for(int d = 0; d < (int)to_delete.size(); d++) {
					map->data[MAP_INDEX(to_delete[d].i, to_delete[d].j, to_delete[d].k)] = level - 1;
				}
				return true;
			}

			printf("Keeping blob with size = %d voxels\n", num_filled);
			return false;
		}

		bool flood_fill_cavity(int i, int j, int k, map_data *map, float level, int *workspace, int fill_threshold) {
			std::vector<map_pos> to_check;
			std::vector<map_pos> to_fill;
			to_check.push_back(map_pos(i,j,k));
			int num_filled = 0;
			while(to_check.size() > 0) {
				map_pos p = to_check.back();
				if(workspace[MAP_INDEX(p.i, p.j, p.k)] != 0) {
					to_check.pop_back();
				} else {
					if(map->data[MAP_INDEX(p.i, p.j, p.k)] <= level) {
						workspace[MAP_INDEX(p.i, p.j, p.k)] = 1;
						to_fill.push_back(map_pos(p.i, p.j, p.k));
						num_filled++;
						to_check.push_back(map_pos(p.i-1,p.j,p.k));
						to_check.push_back(map_pos(p.i+1,p.j,p.k));
						to_check.push_back(map_pos(p.i,p.j-1,p.k));
						to_check.push_back(map_pos(p.i,p.j+1,p.k));
						to_check.push_back(map_pos(p.i,p.j,p.k-1));
						to_check.push_back(map_pos(p.i,p.j,p.k+1));
					} else {
						to_check.pop_back();
					}
				}
			}

			if(num_filled <= fill_threshold) {
				for(int d = 0; d < (int)to_fill.size(); d++) {
					map->data[MAP_INDEX(to_fill[d].i, to_fill[d].j, to_fill[d].k)] = level + 1;
				}
				return true;
			}

			printf("Keeping cavity with size = %d voxels\n", num_filled);
			return false;
		}
};

// Extract the header and voxel data from the given emdb binary map file
int extract_data(char *map_fname, header_info *header, map_data *map)
{
	FILE *mapfile = NULL;
	int i;
	int shutup;

	if((mapfile = fopen(map_fname, "r")) == NULL) {
		printf("Could not open map file %s\n", map_fname);
		return 1;
	}

	shutup = fread(&(header)->NC, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NR, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NS, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->MODE, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NCSTART, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NRSTART, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NSSTART, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NX, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NY, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NZ, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->X_LENGTH, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->Y_LENGTH, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->Z_LENGTH, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->ALPHA,  sizeof(float), 1, mapfile);
	shutup = fread(&(header)->BETA, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->GAMMA, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->MAPC, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->MAPR, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->MAPS, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->AMIN, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->AMAX, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->AMEAN, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->ISPG, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NSYMBT, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->LSKFLG, sizeof(int32_t), 1, mapfile);

	for(i = 0; i < 9; i++) {
		shutup = fread(&(header->SKWMAT[i]), sizeof(float), 1, mapfile);
	}

	for(i = 0; i < 3; i++) {
		shutup = fread(&(header->SKWTRN[i]), sizeof(float), 1, mapfile);
	}

	for(i = 0; i < 15; i++) {
		shutup = fread(&(header->EXTRA[i]), sizeof(int32_t), 1, mapfile);
	}

	header->MAP[0] = 'B'; header->MAP[1] = 'A'; header->MAP[2] = 'D'; header->MAP[3] = '\0';
	shutup = fread(header->MAP, sizeof(char), 4, mapfile);
	header->MAP[4] = '\0';

	shutup = fread(&(header->MACHST), sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header->RMS), sizeof(float), 1, mapfile);
	shutup = fread(&(header->NLABL), sizeof(int32_t), 1, mapfile);

	shutup = fread(header->LABEL_N, sizeof(char), 200 * WORD_SIZE, mapfile);
	header->LABEL_N[200 * WORD_SIZE] = '\0';

//	shutup = fread(header->LABEL_N, sizeof(char), 1, mapfile);

	// Now read in the density data
/*
	map->num_voxels = header->NX * header->NY * header->NZ;
	if((header->NX != header->NY) || (header->NX != header->NZ)) {
		printf("Error: NX, NY and NZ should all be the same\n");
		return 1;
	}
	map->N = header->NX;
	map->voxel_size = header->X_LENGTH/(float)map->N;
	map->data = new float[map->num_voxels];
	shutup = fread(map->data, sizeof(float), map->num_voxels, mapfile);
*/

	if((header->NX < 1) || (header->NY < 1) || (header->NZ < 1)) {
		printf("Error: Map file has dimensions %dx%dx%d\n", header->NX, header->NY, header->NZ);
		return 1;
	}

	map->num_voxels = (header->NX + 2) * (header->NY + 2) * (header->NZ + 2);
//	if((header->NX != header->NY) || (header->NX != header->NZ)) {
//		printf("Error: NX, NY and NZ should all be the same\n");
//		return 1;
//	}
	map->NX = header->NX + 2;
	map->NY = header->NY + 2;
	map->NZ = header->NZ + 2;
	map->voxel_sizeX = header->X_LENGTH/(float)(map->NX - 2);
	map->voxel_sizeY = header->Y_LENGTH/(float)(map->NY - 2);
	map->voxel_sizeZ = header->Z_LENGTH/(float)(map->NZ - 2);
	map->data = new float[map->num_voxels];
	for(i = 0; i < map->num_voxels; i++) {
		map->data[i] = -INFINITY;
	}
	for(int z = 1; z < map->NZ - 1; z++) {
		for(int y = 1; y < map->NY - 1; y++) {
			for(i = MAP_INDEX(1,y,z); i < MAP_INDEX(header->NX + 1,y,z); i++) {
				shutup = fread(&map->data[i], sizeof(float), 1, mapfile);
			}
		}
	}

	fclose(mapfile);

	return 0;
}

void calc_center_position(char *f, double *x, double *y, double *z) {
	FILE *fin;	
	if((fin = fopen(f, "r")) == NULL) {
		printf("Could not open pdb file '%s'\n", f);
		return;
	}
	const int max_line_length = 100;
	char line[max_line_length];
	char record_name[5];
	int i = 0;
	double xmax, xmin, ymax, ymin, zmax, zmin, temp;
	xmin = INFINITY;
	xmax = -INFINITY;
	ymin = INFINITY;
	ymax = -INFINITY;
	zmin = INFINITY;
	zmax = -INFINITY;
	*x = 0.0;
	*y = 0.0;
	*z = 0.0;
	for(;;) {
		if(fgets(line, max_line_length, fin) != NULL) {

			// get first 4 characters of line
			for(int j = 0; j < 5; j++) {
				record_name[j] = line[j];
			}
			if(strstr(record_name, "ATOM ") != NULL) {
				char sx[9], sy[9], sz[9];
				int check = sscanf(line, "%*30c%8c%8c%8c", sx, sy, sz);
				if(check != 3) {
					Error_text();
					printf("Couldn't read line in frame at atom %d\n", i);
					*x = 0.0;
					*y = 0.0;
					*z = 0.0;
				}

				// work around for faulty 8 width field handling
				sx[8] = '\0';
				sy[8] = '\0';
				sz[8] = '\0';
				temp = atof(sx);
				if(temp < xmin) {
					xmin = temp;		
				} 
				if(temp > xmax) {
					xmax = temp;
				}

				temp = atof(sy);
				if(temp < ymin) {
					ymin = temp;		
				} 
				if(temp > ymax) {
					ymax = temp;
				}

				temp = atof(sz);
				if(temp < zmin) {
					zmin = temp;		
				} 
				if(temp > zmax) {
					zmax = temp;
				}
				i++;
			} 
		} else {
			printf("End of .pdb reached. %d atoms read\n", i);
			*x = (xmax + xmin)/2.0;
			*y = (ymax + ymin)/2.0;
			*z = (zmax + zmin)/2.0;

			return;
		}
	}
}
void print_data_to_file(FILE *f, map_data *map) {
	int i;
	fprintf(f, "#\n# Voxel density data:\n");
	for(i = 0; i < map->num_voxels; i++) {
		fprintf(f, "%f\n", map->data[i]);
	}
}

void print_filtered_data_to_file(FILE *f, map_data *map, float level) {
	int i;
	fprintf(f, "#\n# Voxels filtered at level %f:\n", level);
	for(i = 0; i < map->num_voxels; i++) {
		if(map->data[i] < level) {
			fprintf(f, "0\n");
		} else {
			fprintf(f, "1\n");
		}
	}
}

void print_header_to_file(FILE *f, header_info *header)
{
	int i;

	fprintf(f, "#\n# Header:\n");
	fprintf(f, "# NC = %d\n", header->NC);
	fprintf(f, "# NR = %d\n", header->NR);
	fprintf(f, "# NS = %d\n", header->NS);
	fprintf(f, "# MODE = %d\n", header->MODE);
	fprintf(f, "# NCSTART = %d\n", header->NCSTART);
	fprintf(f, "# NRSTART = %d\n", header->NRSTART);
	fprintf(f, "# NSSTART = %d\n", header->NSSTART);
	fprintf(f, "# NX = %d\n", header->NX);
	fprintf(f, "# NY = %d\n", header->NY);
	fprintf(f, "# NZ = %d\n", header->NZ);
	fprintf(f, "# X_LENGTH = %f\n", header->X_LENGTH);
	fprintf(f, "# Y_LENGTH = %f\n", header->Y_LENGTH);
	fprintf(f, "# Z_LENGTH = %f\n", header->Z_LENGTH);
	fprintf(f, "# ALPHA = %f\n", header->ALPHA);
	fprintf(f, "# BETA = %f\n", header->BETA);
	fprintf(f, "# GAMMA = %f\n", header->GAMMA);
	fprintf(f, "# MAPC = %d\n", header->MAPC);
	fprintf(f, "# MAPR = %d\n", header->MAPR);
	fprintf(f, "# MAPS = %d\n", header->MAPS);
	fprintf(f, "# AMIN = %f\n", header->AMIN);
	fprintf(f, "# AMAX = %e\n", header->AMAX);
	fprintf(f, "# AMEAN = %e\n", header->AMEAN);
	fprintf(f, "# ISPG = %d\n", header->ISPG);
	fprintf(f, "# NSYMBT = %d\n", header->NSYMBT);
	fprintf(f, "# LSKFLG = %d\n", header->LSKFLG);
	fprintf(f, "# SKWMAT =\n# ");
	for(i = 0; i < 9; i++) {
		fprintf(f, "%f ", header->SKWMAT[i]);
		if((i+1)%3 == 0) {
			fprintf(f, "\n# ");
		}
	}
	fprintf(f, "# SKWTRN =\n# ");
	for(i = 0; i < 3; i++) {
		fprintf(f, "%f ", header->SKWTRN[i]);
	}
	fprintf(f, "\n");
	fprintf(f, "# EXTRA =\n# ");
	for(i = 0; i < 14; i++) {
		fprintf(f, "%d ", header->EXTRA[i]);
	}
	fprintf(f, "\n");
	fprintf(f, "# MAP = %s\n", header->MAP);
	fprintf(f, "# MACHST = %d\n", header->MACHST);
	fprintf(f, "# RMS = %f\n", header->RMS);
	fprintf(f, "# NLABL = %d\n", header->NLABL);
	fprintf(f, "# LABEL_N = %s\n", header->LABEL_N);
}

void coarsen_map(map_data *map, int coarseness)
{
	// Find the new map dimensions after coarsenning
	int new_NX = (int)floor((map->NX - 2) / (float)coarseness) + 2;
	int new_NY = (int)floor((map->NY - 2) / (float)coarseness) + 2;
	int new_NZ = (int)floor((map->NZ - 2) / (float)coarseness) + 2;
	int new_num_voxels = new_NX * new_NY * new_NZ;

	fprintf(stderr, "New dimensions: new_NX = %d, new_NY = %d, new_NZ = %d\n", new_NX, new_NY, new_NZ);	
	// Allocate the new map some memory
	float *new_map = new float[new_num_voxels];
	for(int i = 0; i < new_num_voxels; i++) {
		new_map[i] = -INFINITY;
	}

	float density_sum;
	int num_voxels_per_group = pow(coarseness, 3);
	int index_x, index_y, index_z;
	for(int z = 1; z < new_NZ - 1; z++) {
		for(int y = 1; y < new_NY - 1; y++) {
			for(int x = 1; x < new_NX - 1; x++) {
				density_sum = 0;
				for(int vz = 0; vz < coarseness; vz++) {
					for(int vy = 0; vy < coarseness; vy++) {
						for(int vx = 0; vx < coarseness; vx++) {
							index_x = (x - 1) * coarseness + vx + 1;
							index_y = (y - 1) * coarseness + vy + 1;
							index_z = (z - 1) * coarseness + vz + 1;
							if(index_x < map->NX - 1 && index_y < map->NY - 1 && index_z < map->NZ - 1) {
								density_sum += map->data[MAP_INDEX(index_x, index_y, index_z)];
								//fprintf(stderr, "%e ", density_sum);
							}
						}
					}
				}
				new_map[z * new_NY * new_NX + y * new_NX + x] = density_sum / num_voxels_per_group;
			}
		}
	}

	// Delete the old map data and replace it with the new, coarser version
	map->NX = new_NX;
	map->NY = new_NY;
	map->NZ = new_NZ;
	map->num_voxels = new_num_voxels;
	map->voxel_sizeX *= coarseness;
	map->voxel_sizeY *= coarseness;
	map->voxel_sizeZ *= coarseness;
	delete[] map->data;
	map->data = new_map;
}

void coarsen_map_2(map_data *map, int coarseness)
{
	// Find the new map dimensions after coarsenning
	int coarse_x = coarseness, coarse_y = coarseness, coarse_z = coarseness;
	int new_NX;
	int new_NY;
	int new_NZ;
	while((map->NX - 2) % coarse_x != 0) {
		coarse_x++;
	}
	while((map->NY - 2) % coarse_y != 0) {
		coarse_y++;
	}
	while((map->NZ - 2) % coarse_z != 0) {
		coarse_z++;
	}
	new_NX = ((map->NX - 2) / coarse_x) + 2;
	new_NY = ((map->NY - 2) / coarse_y) + 2;
	new_NZ = ((map->NZ - 2) / coarse_z) + 2;
	int new_num_voxels = new_NX * new_NY * new_NZ;
	fprintf(stderr, "Rounded coarseness: coarse_x = %d, coarse_y = %d, coarse_z = %d\n", coarse_x, coarse_y, coarse_z);
	fprintf(stderr, "New dimensions: new_NX = %d, new_NY = %d, new_NZ = %d\n", new_NX, new_NY, new_NZ);	
	// Allocate the new map some memory
	float *new_map = new float[new_num_voxels];
	for(int i = 0; i < new_num_voxels; i++) {
		new_map[i] = -INFINITY;
	}

	float density_sum;
	int num_voxels_per_group = coarse_x * coarse_y * coarse_z;

	for(int z = 1; z < new_NZ - 1; z++) {
		for(int y = 1; y < new_NY - 1; y++) {
			for(int x = 1; x < new_NX - 1; x++) {
				density_sum = 0;
				for(int vz = 0; vz < coarse_z; vz++) {
					for(int vy = 0; vy < coarse_y; vy++) {
						for(int vx = 0; vx < coarse_x; vx++) {
							density_sum += map->data[MAP_INDEX((x - 1) * coarse_x + vx + 1, (y - 1) * coarse_y + vy + 1, (z - 1) * coarse_z + vz + 1)];
							//fprintf(stderr, "%e ", density_sum);
						}
					}
				}
				new_map[z * new_NY * new_NX + y * new_NX + x] = density_sum / num_voxels_per_group;
			}
		}
	}

	// Delete the old map data and replace it with the new, coarser version
	map->NX = new_NX;
	map->NY = new_NY;
	map->NZ = new_NZ;
	map->num_voxels = new_num_voxels;
	map->voxel_sizeX *= coarse_x;
	map->voxel_sizeY *= coarse_y;
	map->voxel_sizeZ *= coarse_z;
	delete[] map->data;
	map->data = new_map;

}

int get_euler_number(map_data *map, float level)
{
	int chi = 0;
	int flags[8];
	for(int z = 0; z < map->NZ - 1; z++) {
		for(int y = 0; y < map->NY - 1; y++) {
			for(int x = 0; x < map->NX - 1; x++) {
				flags[0] = (map->data[MAP_INDEX(x    , y,     z    )] <= level);
				flags[1] = (map->data[MAP_INDEX(x + 1, y,     z    )] <= level);
				flags[2] = (map->data[MAP_INDEX(x + 1, y + 1, z    )] <= level);
				flags[3] = (map->data[MAP_INDEX(x    , y + 1, z    )] <= level);
				flags[4] = (map->data[MAP_INDEX(x,     y,     z + 1)] <= level);
				flags[5] = (map->data[MAP_INDEX(x + 1, y,     z + 1)] <= level);
				flags[6] = (map->data[MAP_INDEX(x + 1, y + 1, z + 1)] <= level);
				flags[7] = (map->data[MAP_INDEX(x    , y + 1, z + 1)] <= level);
				chi += get_euler_characteristic_x_4(flags);
			}
		}
	}
	return chi/4;
}

int write_map_to_file(map_data *map, header_info *header, char *fname)
{
	FILE *out = NULL;
	if((out = fopen(fname, "w")) == NULL) {
		printf("Could not open %s for writing.\n", fname);
		return -1;
	}

	// write the header to file
	fwrite(&(header)->NC, sizeof(int32_t), 1, out);
	fwrite(&(header)->NR, sizeof(int32_t), 1, out);
	fwrite(&(header)->NS, sizeof(int32_t), 1, out);
	fwrite(&(header)->MODE, sizeof(int32_t), 1, out);
	fwrite(&(header)->NCSTART, sizeof(int32_t), 1, out);
	fwrite(&(header)->NRSTART, sizeof(int32_t), 1, out);
	fwrite(&(header)->NSSTART, sizeof(int32_t), 1, out);
	fwrite(&(header)->NX, sizeof(int32_t), 1, out);
	fwrite(&(header)->NY, sizeof(int32_t), 1, out);
	fwrite(&(header)->NZ, sizeof(int32_t), 1, out);
	fwrite(&(header)->X_LENGTH, sizeof(float), 1, out);
	fwrite(&(header)->Y_LENGTH, sizeof(float), 1, out);
	fwrite(&(header)->Z_LENGTH, sizeof(float), 1, out);
	fwrite(&(header)->ALPHA,  sizeof(float), 1, out);
	fwrite(&(header)->BETA, sizeof(float), 1, out);
	fwrite(&(header)->GAMMA, sizeof(float), 1, out);
	fwrite(&(header)->MAPC, sizeof(int32_t), 1, out);
	fwrite(&(header)->MAPR, sizeof(int32_t), 1, out);
	fwrite(&(header)->MAPS, sizeof(int32_t), 1, out);
	fwrite(&(header)->AMIN, sizeof(float), 1, out);
	fwrite(&(header)->AMAX, sizeof(float), 1, out);
	fwrite(&(header)->AMEAN, sizeof(float), 1, out);
	fwrite(&(header)->ISPG, sizeof(int32_t), 1, out);
	fwrite(&(header)->NSYMBT, sizeof(int32_t), 1, out);
	fwrite(&(header)->LSKFLG, sizeof(int32_t), 1, out);

	for(int i = 0; i < 9; i++) {
		fwrite(&(header->SKWMAT[i]), sizeof(float), 1, out);
	}

	for(int i = 0; i < 3; i++) {
		fwrite(&(header->SKWTRN[i]), sizeof(float), 1, out);
	}

	for(int i = 0; i < 15; i++) {
		fwrite(&(header->EXTRA[i]), sizeof(int32_t), 1, out);
	}
	fwrite(header->MAP, sizeof(char), 4, out);

	fwrite(&(header->MACHST), sizeof(int32_t), 1, out);
	fwrite(&(header->RMS), sizeof(float), 1, out);
	fwrite(&(header->NLABL), sizeof(int32_t), 1, out);

	fwrite(header->LABEL_N, sizeof(char), 200 * WORD_SIZE, out);

	for(int z = 1; z < map->NZ - 1; z++) {
		for(int y = 1; y < map->NY - 1; y++) {
			for(int x = 1; x < map->NX - 1; x++) {
				fwrite(&(map->data[MAP_INDEX(x, y, z)]), sizeof(float), 1, out);
			}
		}
	}

	return 0;
}

int main(int argc, char **argv)
{
	int arg;
	char map_fname[MAX_FNAME_LENGTH], out_fname[MAX_FNAME_LENGTH], out_format[MAX_FNAME_LENGTH], pdb_fname[MAX_FNAME_LENGTH];
	float level = 0;
	int coarseness = 0;
	int choose_format = OUTPUT_FORMAT_TEXT;
	int interpolate = 1;
	int cull_threshold = -1;
	int fill_threshold = -1;

	int set_map_fname = 0, set_out_fname = 0, set_level = 0, set_coarseness = 0, set_format = 0, set_euler = 0, set_pdb_fname = 0, set_cull_floaters = 0, set_fill_cavities = 0;
	
	double cent_x, cent_y, cent_z;

	if(argc == 1) {
		printf("Usage: ./emdb_map_to_ffea -map INPUT_MAP_FNAME -out OUTPUT_FNAME -pdb ORIGINAL PDB [-level LEVEL FILTER] [-format OUTPUT FORMAT] [-coarse COARSENING LEVEL] [-interpolate yes/no] [-euler] [-cull_floaters THRESHOLD IN VOXELS] [-fill_cavities [THRESHOLD IN VOXELS]]\n");
		return 1;
	}

	arg = 1;
	while(arg < argc) {
		if(strcmp(argv[arg], "-map") == 0) {
			arg++;
			if(arg < argc) {
				sprintf(map_fname, "%s", argv[arg]);
				printf("Input map file name = %s\n", map_fname);
				set_map_fname = 1;
			} else {
				Error_text();
				printf("option -map is missing argument (Usage: -map INPUT_MAP_FILENAME)\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-out") == 0) {
			arg++;
			if(arg < argc) {
				sprintf(out_fname, "%s", argv[arg]);
				printf("Output text file name = %s\n", out_fname);
				set_out_fname = 1;
			} else {
				Error_text();
				printf("option -out is missing argument (Usage: -out OUTPUT_TEXT_FILENAME)\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-level") == 0) {
			arg++;
			if(arg < argc) {
				level = atof(argv[arg]);
				printf("Filter level = %f\n", level);
				set_level = 1;
			} else {
				Error_text();
				printf("option -level is missing argument (Usage: -level FILTER_LEVEL)\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-coarse") == 0) {
			arg++;
			if(arg < argc) {
				coarseness = atoi(argv[arg]);
				printf("Coarseness = %d\n", coarseness);
				set_coarseness = 1;
			} else {
				Error_text();
				printf("option -coarse is missing argument (Usage: -coarse COARSENING_FACTOR)\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-interpolate") == 0) {
			arg++;
			if(arg < argc) {
				if(strcmp(argv[arg], "yes") == 0) {
					interpolate = 1;
				} else {
					interpolate = 0;
				}
				printf("interpolate = %d\n", interpolate);
			} else {
				Error_text();
				printf("option -interpolate is missing argument (Usage: -interpolate 'yes' or 'no')\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-format") == 0) {
			arg++;
			if(arg < argc) {
				sprintf(out_format, "%s", argv[arg]);

				// check if format is recognised
				if(strcmp(out_format, "text") == 0) {
					printf("Output format = text file\n");
					choose_format = OUTPUT_FORMAT_TEXT;
				} else if(strcmp(out_format, "surf") == 0) {
					printf("Output format = netgen surf\n");
					choose_format = OUTPUT_FORMAT_SURF;
				} else if(strcmp(out_format, "obj") == 0) {
					printf("Output format = wavefront obj\n");
					choose_format = OUTPUT_FORMAT_OBJ;
				} else if(strcmp(out_format, "stl") == 0) {
					printf("Output format = stereolithography stl\n");
					choose_format = OUTPUT_FORMAT_STL;
				} else if(strcmp(out_format, "off") == 0) {
					printf("Output format = Geomview Object\n");
					choose_format = OUTPUT_FORMAT_OFF;
				} else if(strcmp(out_format, "map") == 0) {
					printf("Output format = CCP4 density map\n");
					choose_format = OUTPUT_FORMAT_MAP;
				} else {
					Error_text();
					printf("Unrecognised format '%s'. Recognised formats are:\nsurf\nobj\nstl\noff\nmap\n", out_format);
					return 1;
				}

				set_format = 1;
			} else {
				Error_text();
				printf("option -format is missing argument (Usage: -format OUTPUT_FORMAT)\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-cull_floaters") == 0) {
			arg++;
			if(arg < argc) {
				set_cull_floaters = 1;
				cull_threshold = atoi(argv[arg]);
				printf("cull_threshold = %d\n", cull_threshold);
			} else {
				Error_text();
				printf("option -cull_floaters is missing argument (Usage: -cull_floaters [THRESHOLD IN NUMBER OF VOXELS])\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-fill_cavities") == 0) {
			arg++;
			if(arg < argc) {
				set_fill_cavities = 1;
				fill_threshold = atoi(argv[arg]);
				printf("fill_threshold = %d\n", fill_threshold);
			} else {
				Error_text();
				printf("option -fill_cavities is missing argument (Usage: -fill_cavities [THRESHOLD IN NUMBER OF VOXELS])\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-pdb") == 0) {
			arg++;
			if(arg < argc) {
				set_pdb_fname = 1;
				sprintf(pdb_fname, "%s", argv[arg]);
				printf("pdb_fname = %s\n", pdb_fname);
			} else {
				Error_text();
				printf("option -pdb is missing argument (Usage: -pdb [PDB FILENAME MAP WAS CREATED FROM])\n");
				return 1;
			}
		} else if(strcmp(argv[arg], "-euler") == 0) {
			set_euler = 1;
			printf("calculate euler characteristic number -> yes\n");
		} else if(strcmp(argv[arg], "help") == 0) {
			printf("./emdb_map_to_ffea\n");
			printf("Recognised options:\n");
			printf("-map [Followed by the name of the input binary map file]\n");
			printf("-out [Followed by the desired name of the output file]\n");
			printf("-level [Followed by the density level at which to filter the map]\n");
			printf("-format [Followed by: text, surf, obj, stl or off]\n");
			printf("-coarse [Followed by an integer, N, causing groups of NxNxN voxels to be averaged]\n");
			printf("-interpolate [Followed by yes or no. If yes, then linearly interpolate point position based on density]\n");
			printf("-pdb [Followed by a pdb filename, inc. extention, if map was created from a pdb file]\n");
			printf("-euler [In order to calculate the euler characteristic number of the map (after any coarsening)]\n");
		} else {
			Error_text();
			printf("Unrecgonised option '%s'\n", argv[arg]);
			printf("For info use 'emdb_map_to_ffea help'\n");
			return 1;
		}

		arg++;
	}

	// check that user has provided input and output file names
	if(set_map_fname == 0) {
		Error_text();
		printf("input map file not specified. (Set this using '-map' option)\n");
		return 1;
	}

	if(set_out_fname == 0) {
		Error_text();
		printf("No form of output specified. (Use '-out' to produce a text file or a surface mesh)\n");
		return 1;
	}

	if(set_format == 0) {
		Error_text();
		printf("Output format not specified. (Set this using '-format' option)\n");
		return 1;
	}

	if(choose_format == OUTPUT_FORMAT_SURF || choose_format == OUTPUT_FORMAT_OFF || choose_format == OUTPUT_FORMAT_OBJ || choose_format == OUTPUT_FORMAT_STL) {
		if(set_level == 0) {
			Error_text();
			printf("Cannot create surface mesh with density filter isolevel unset (use '-level' to set it)\n");
			return 1;
		}
	}

	// extract header and voxel data from input file
	header_info header;
	map_data map;
	if(extract_data(map_fname, &header, &map) != 0) {
		Error_text();
		printf("Couldn't extract data from map file.\n");
		return 1;
	}

	// Print so user knows what the binary data actually looks like...
	printf("%s header:\n", map_fname);
	print_header_to_file(stdout, &header);

	// If selected, coarsen the mesh to the specified degree
	if(set_coarseness == 1) {
		printf("Coarsenning map (coarseness level %d)\n", coarseness);
		coarsen_map(&map, coarseness);
	}

	// fill any cavities with sizes below the given threshold
	if(set_fill_cavities == 1) {
		map_flood_filler mff;
		mff.fill_cavities(&map, level, fill_threshold);
	}

	// Cull any 'floaters' (below threshold number of voxels)
	if(set_cull_floaters == 1) {
		map_flood_filler mff;
		mff.cull_floaters(&map, level, cull_threshold);
	}

	// If chosen, will use pdb file to recalculate center that map creation lost!
	if(set_pdb_fname == 1) {
		calc_center_position(pdb_fname, &cent_x, &cent_y, &cent_z);
		printf("Centroid = (%f %f %f)\n", cent_x, cent_y, cent_z);
		cent_x = cent_x - header.X_LENGTH/2.0;
		cent_y = cent_y - header.Y_LENGTH/2.0;
		cent_z = cent_z - header.Z_LENGTH/2.0;
		printf("->Shift = (%f %f %f)\n", cent_x, cent_y, cent_z);
	} else {
		printf("Assuming no pdb file. No centroid will be calculated\n");
	}

	// If MAP format selected, output the map (after all the coarsening, filling etc.)
	if(choose_format == OUTPUT_FORMAT_MAP) {
		printf("Writing CCP4 map file...\n");
		header.NC = map.NX - 2;
		header.NR = map.NY - 2;
		header.NS = map.NZ - 2;
		header.NX = map.NX - 2;
		header.NY = map.NY - 2;
		header.NZ = map.NZ - 2;
		header.X_LENGTH = header.NX * map.voxel_sizeX;
		header.Y_LENGTH = header.NY * map.voxel_sizeY;
		header.Z_LENGTH = header.NZ * map.voxel_sizeZ;
		write_map_to_file(&map, &header, out_fname);
	}

	// print this to a text file
	if(choose_format == OUTPUT_FORMAT_TEXT) {
		FILE *outfile = NULL;
		if((outfile = fopen(out_fname, "w")) == NULL) {
			printf("Could not open output file %s for writing\n", out_fname);
			return 1;
		}
		fprintf(outfile, "# emdb_map_to_text output\n#\n# Options:\n# -map %s\n# -out %s\n", map_fname, out_fname);
		if(set_level == 1) {
			fprintf(outfile, "# -level %f\n", level);
		}
		print_header_to_file(outfile, &header);

		if(set_level == 0) {
			print_data_to_file(outfile, &map);
		} else {
			print_filtered_data_to_file(outfile, &map, level);
		}
		fprintf(outfile, "\n");
		fclose(outfile);
		return 0;
	}

	// Generate an output surface file if requested
	if(choose_format == OUTPUT_FORMAT_SURF || choose_format == OUTPUT_FORMAT_OFF || choose_format == OUTPUT_FORMAT_OBJ || choose_format == OUTPUT_FORMAT_STL) {
		map_to_surf_generator m2sgen;
		printf("Using marching cubes to generate map surface...\n");
		m2sgen.generate_surface(&map, level, interpolate);

		if(choose_format == OUTPUT_FORMAT_SURF) {
			printf("Writing netgen surf file...\n");
			m2sgen.write_netgen_surf_file(out_fname, cent_x, cent_y, cent_z);
		} else if(choose_format == OUTPUT_FORMAT_OBJ) {
			printf("Writing OBJ file...\n");
			m2sgen.write_OBJ_file(out_fname, cent_x, cent_y, cent_z);
		} else if(choose_format == OUTPUT_FORMAT_STL) {
			printf("Writing STL file...\n");
			m2sgen.write_STL_file(out_fname, cent_x, cent_y, cent_z);
		} else if(choose_format == OUTPUT_FORMAT_OFF) {
			printf("Writing OFF file...\n");
			m2sgen.write_OFF_file(out_fname, cent_x, cent_y, cent_z);
		}
	}


	// if chosen, calc the euler characteristic number of the density map
	if(set_euler == 1) {
		printf("Calculating euler characteristic number of map...\n");
		printf("Euler characteristic = %d\n", get_euler_number(&map, level));
	}

	return 0;
}
