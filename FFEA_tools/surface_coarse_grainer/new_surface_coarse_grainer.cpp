#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <string.h>

class vector3
{
	public:
		vector3() {
			x = 0.0;
			y = 0.0;
			z = 0.0;
		}
		
		vector3(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
		}

		void normalise()
		{
			double l = sqrt(x * x + y * y + z * z);
			x /= l;
			y /= l;
			z /= l;
		}

		void midpoint(vector3 a, vector3 b) 
		{
			this->x = (a.x + b.x)/2.0;
			this->y = (a.y + b.y)/2.0;
			this->z = (a.z + b.z)/2.0;	
		}

		double x, y, z;
};

class Face
{
	public:
		Face() {
			this->n1 = 0;
			this->n2 = 0;
			this->n3 = 0;
			this->remaining_node = -100;
		}

		Face(int n1, int n2, int n3) {
			this->n1 = n1;
			this->n2 = n2;
			this->n3 = n3;
			this->remaining_node = -100;
		}

		Face(Face *f) {
			this->n1 = f->n1;
			this->n2 = f->n2;
			this->n3 = f->n3;
			this->remaining_node = f->remaining_node;
		}

		void calc_shortest_side(vector3 *nodes)
		{
			// Get the lengths of each of the face's edges
			double	l12 = edge_length(&nodes[n1], &nodes[n2]),
				l13 = edge_length(&nodes[n1], &nodes[n3]),
				l23 = edge_length(&nodes[n2], &nodes[n3]);

			// Work out which edge is the shortest
			if(l12 < l13) {
				if(l12 < l23) {
					// l12 is shortest
					shortest_length = l12;
					shortest_edge_node_a = n1;
					shortest_edge_node_b = n2;
					remaining_node = n3;
				} else {
					// l23 is shortest
					shortest_length = l23;
					shortest_edge_node_a = n2;
					shortest_edge_node_b = n3;
					remaining_node = n1;
				}
			} else {
				if(l13 < l23) {
					// l13 is shortest
					shortest_length = l13;
					shortest_edge_node_a = n1;
					shortest_edge_node_b = n3;
					remaining_node = n2;
				} else {
					// l23 is shortest
					shortest_length = l23;
					shortest_edge_node_a = n2;
					shortest_edge_node_b = n3;
					remaining_node = n1;
				}
			}
		}

		int get_node(int n)
		{
			if(n == 1) {
				return n1;
			} else if(n == 2) {
				return n2;
			} else if(n == 3) {
				return n3;
			} else {
				return -1;
			}
		}

		int get_shortest_node(int n) 
		{
			if(n == 1) {
				return shortest_edge_node_a;
			} else if (n == 2) {
				return shortest_edge_node_b;
			} else {
				return -1;
			}
		}
		
		int get_remaining_node()
		{
			return remaining_node;
		}

		double get_shortest_length() 
		{
			return shortest_length;	
		}

		bool do_i_have_this_edge(int node_a, int node_b)
		{	
			if(n1 == node_a) {
				if(n2 == node_b || n3 == node_b) {
					return true;
				}
			} 

			if(n2 == node_a) {
				if(n1 == node_b || n3 == node_b) {
					return true;
				}
			} 
			
			if(n3 == node_a) {
				if(n1 == node_b || n2 == node_b) {
					return true;
				}
			}

			return false;
		}

		bool do_i_have_this_node(int node_a) 
		{
			if(n1 == node_a || n2 == node_a || n3 == node_a) {
				return true;
			} else {
				return false;
			}
		}

		void reconnect_from_to(int from, int to) {

			if(n1 == from) {
				n1 = to;
				return;
			} else if(n2 == from) {
				n2 = to;
				return;
			} else if(n3 == from) {
				n3 = to;
				return;
			} 
		}

		void remap_node_indices(int *map)
		{
			n1 = map[n1];
			n2 = map[n2];
			n3 = map[n3];
		}

		private:
			int n1, n2, n3;
			int remaining_node;
			double shortest_length;
			int shortest_edge_node_a, shortest_edge_node_b;

			double edge_length(vector3 *a, vector3 *b)
			{
				double x_length = a->x - b->x, y_length = a->y - b->y, z_length = a->z - b->z;
				return sqrt(pow(x_length, 2) + pow(y_length, 2) + pow(z_length, 2));
			}
};

class Surface 
{
	public:
		int init(char *insurf_fname)
		{
			printf("Reading in surface file %s...\n", insurf_fname);
			// Open surface file
			FILE *infile;
			if((infile = fopen(insurf_fname, "r")) == NULL) {
				printf("Error while trying to open %s for reading. Does file exist?\n", insurf_fname);
				return -1;
			}

			// Check that file is a netgen surface file
			char surfacemesh[13];
			if(fgets(surfacemesh, 12, infile) == NULL) {
				printf("Error when reading first line of %s.\n", insurf_fname);
				return -1;
			}
			if(strcmp(surfacemesh, "surfacemesh") != 0) {
				printf("Error, not a netgen surf file. First line of %s should read 'surfacemesh', not '%s'.\n", insurf_fname, surfacemesh);
				return -1;
			}

			// get the number of nodes
			if(fscanf(infile, "%d\n", &num_nodes) != 1) {
				printf("Error when reading number of nodes line of %s.\n", insurf_fname);
				return -1;
			}

			// remember how many nodes there were originally
			num_nodes_original = num_nodes;

			// read nodes into array
			nodes = new vector3[num_nodes];
			for(int i = 0; i < num_nodes; i++) {
				if(fscanf(infile, "%lf %lf %lf", &nodes[i].x, &nodes[i].y, &nodes[i].z) != 3) {
					printf("Error when reading node %d.\n", i);
					return -1;
				}
			}

			// allocate memory for the node index mapping (calculated at the end of coarsen())
			node_map = new int[num_nodes];

			// get the number of faces
			if(fscanf(infile, "%d\n", &num_faces) != 1) {
				printf("Error when reading number of faces line of %s.\n", insurf_fname);
				return -1;
			}

			// remember how many faces there were originally
			num_faces_original = num_faces;

			// read faces into linked list
			int n1, n2, n3, n4, n5, n6;
			for(int i = 0; i < num_faces; i++) {
				if(fscanf(infile, "%d %d %d", &n1, &n2, &n3) != 3 ) {
					printf("Error when reading face %d.\n %d %d %d\n", i, n1, n2, n3);
					return -1;
				}
				faces.push_back(new Face(n1 - 1, n2 - 1, n3 - 1));
			}
			printf("...Done. Read in %d nodes and %d faces.\n", num_nodes, num_faces);

			// calculate the shortest edge of all faces
			printf("Calculating all edge lengths and storing the shortest...");
			double length, min_length = INFINITY;
			for(std::list<Face*>::iterator face_iter = faces.begin() ; face_iter != faces.end(); face_iter++) {
				(*face_iter)->calc_shortest_side(nodes);
				length = (*face_iter)->get_shortest_length();
				if(length < min_length) {
					min_length = length;
				}
			}
			printf("...Done. Current Shortest Edge = %e\n", min_length);
			return 0;
		}

		int coarsen(double thresh)
		{

			std::list<Face*>::iterator face_iter, face_to_delete[4];
			int have_we_finished, success;
			int node_to_delete_a, node_to_delete_b;
			int nodes_deleted = 0, faces_deleted = 0;
			double length, smallest_length;
			vector3 center_position, dead_position(INFINITY, INFINITY, INFINITY);
			int count = 0;
			int node_a, node_b, node_c, test_node;
			while(true) {
				
				have_we_finished = 1;

				// Firstly locate smallest length than the threshold in the list
				smallest_length = thresh;
				for(face_iter = faces.begin() ; face_iter != faces.end(); face_iter++) {
					length = (*face_iter)->get_shortest_length();

					if(length < smallest_length) {
						smallest_length = length;
						face_to_delete[0] = face_iter;
						have_we_finished = 0;
					}
				}
				
				// If edges all greater than thresh, return
				if(have_we_finished == 1) {
					printf("System Coarsened to Length %e\n", thresh);
					break;

				} 

				// Assignment of faces and nodes to delete 
				node_to_delete_a = (*face_to_delete[0])->get_shortest_node(1);
				node_to_delete_b = (*face_to_delete[0])->get_shortest_node(2);
				

				// Progress report
				if(count++ == 20) {
					count = 0;
					printf("Faces left = %d\nNodes left = %d\nDeleting Nodes %d(.x = %e) and %d(.x = %e),  Edge of length %f\n\n", num_faces, num_nodes, node_to_delete_a, nodes[node_to_delete_a].x, node_to_delete_b, nodes[node_to_delete_b].x, smallest_length);
				}

				for(face_iter = faces.begin() ; face_iter != faces.end(); face_iter++) {
					node_a = (*face_iter)->get_node(1);
					node_b = (*face_iter)->get_node(2);
					node_c = (*face_iter)->get_node(3);
					if(nodes[node_a].x == INFINITY || nodes[node_b].x == INFINITY || nodes[node_c].x == INFINITY) {
						printf("I think the structure has become flat. System cannot coarsen further this way. Coarsen to a smaller length than %e. Sorry!\n", smallest_length);
						return -1;
					}
				}

				// Making node_to_delete_b the higher of the two in the list
				int temp;
				if(node_to_delete_b < node_to_delete_a) {
					temp = node_to_delete_b;
					node_to_delete_b = node_to_delete_a;
					node_to_delete_a = temp;
				}

				// Finding second face to delete (two faces share an edge)
				int num_sharing_faces = 1;
				for(face_iter = faces.begin() ; face_iter != faces.end(); face_iter++) {
					if(face_iter == faces_to_delete[0]) {
						continue;
					}
					if((*face_iter)->do_i_have_this_edge(node_to_delete_a, node_to_delete_b) == true) {
						if(num_sharing_faces == 1) {
							face_to_delete[1] = face_iter;
						} else if(num_sharing_faces == 2) {
							face_to_delete[2] = face_iter;
						} else if(num_sharing_faces == 3) {
							face_to_delete[3] = face_iter;
						}
						num_sharing_faces++;
					} 
				}

				if(num_sharing_faces == 1) {
					printf("There has been a problem somewhere. Only one face has the shortest edge\n");
					return -1;
				}

				if(num_sharing_faces == 4) {

					// A tetrahedral protrusion has collapsed, leaving a single face sticking out. Deleting protrusion
					printf("Num sharing faces = %d. Tetrahedral protrusion has collapsed. Deleting\n", num_sharing_faces);
				
					// Finding equivalent faces
					success = 0;
					for(int i = 0; i < 3; ++i) {
						test_node = (*face_to_delete[i])->get_remaining_node();
						for(int j = i + 1; j < 4; ++j) {
							if(test_node == (*face_to_delete[j])->get_remaining_node()) {

								// Deleting the faces
								faces.erase(face_to_delete[j]);	
								faces.erase(face_to_delete[i]);	
								faces_deleted = faces_deleted + 2;
								num_faces = num_faces_original - faces_deleted;
	
								// Flagging extra node as dead
								nodes[test_node] = dead_position;
								nodes_deleted++;
								num_nodes = num_nodes_original - nodes_deleted;

								success = 1;
								break;
							}
						}
						if(success == 1) {
							break;	
						}
					}
					// Next iteration
					continue;
				}
				
				// Delete faces from list so they are not in neighbourhood
				faces.erase(face_to_delete[1]);
				faces.erase(face_to_delete[0]);
				faces_deleted = faces_deleted + 2;
				num_faces = num_faces_original - faces_deleted;

				// Finding all faces connected to nodes_to_delete_b and adding to neighbourhood list
				neighbourhood.clear();
				for(face_iter = faces.begin(); face_iter != faces.end(); face_iter++) {
					if((*face_iter)->do_i_have_this_node(node_to_delete_a) == true || (*face_iter)->do_i_have_this_node(node_to_delete_b) == true) {
						neighbourhood.push_back((*face_iter));
					}
				}

				// What type of neighbourhood do we have? (add later) 
					
				// Make all faces connected to node_to_connect_b instead connect to node_to_delete_a i.e 'deleting node_to_delete_b'
				for(face_iter = neighbourhood.begin(); face_iter != neighbourhood.end(); face_iter++) {
					(*face_iter)->reconnect_from_to(node_to_delete_b, node_to_delete_a);
					//printf("Position of node_a: (%e %e %e)\n", nodes[node_to_delete_a].x, nodes[node_to_delete_a].y, nodes[node_to_delete_a].z);
					//printf("Position of node_b: (%e %e %e)\n", nodes[node_to_delete_b].x, nodes[node_to_delete_b].y, nodes[node_to_delete_b].z);
					
				}

				// Move node_to_delete_a to new position: midpoint of nodes_to_delete + some extra distance if volume is constant (add later)
				center_position.midpoint(nodes[node_to_delete_a], nodes[node_to_delete_b]);
				nodes[node_to_delete_a] = center_position;

				// Move node_to_delete_b to infinity but don't delete; out of system but without having to reassign faces within list each step
				nodes[node_to_delete_b] = dead_position;
				nodes_deleted++;
				num_nodes = num_nodes_original - nodes_deleted;

				// Faces within neighbourhood will need to have smallest edges recalculated
				for(face_iter = neighbourhood.begin(); face_iter != neighbourhood.end(); face_iter++) {
					(*face_iter)->calc_shortest_side(nodes);
					length = (*face_iter)->get_shortest_length();
				}
			}
			printf("All done. Faces deleted = %d, Nodes deleted = %d\n", faces_deleted, nodes_deleted);
			return 0;
		}

		int output(char *outsurf_fname) 
		{
			remap_nodes();
			write_surf_file(outsurf_fname);
			return 0;
		}

		int remap_nodes()
		{
			
			int i, node_count = 0;
			for(i = 0; i < num_nodes_original; ++i) {
				if(nodes[i].x == INFINITY) {
					node_map[i] = -1;
				} else {
					node_map[i] = node_count;
					nodes[node_count] = nodes[i];
					node_count++;
				}
			}
			
			// remap the node indices in the faces array
			for(std::list<Face*>::iterator face_iter = faces.begin(); face_iter != faces.end(); face_iter++) {
				(*face_iter)->remap_node_indices(node_map);
				if((*face_iter)->do_i_have_this_node(-1) == true) {
					printf("Face: (%d %d %d)\n", (*face_iter)->get_node(1), (*face_iter)->get_node(2), (*face_iter)->get_node(3));
					return -1;
				}
			}

			return 0;			
		}

		int write_surf_file(char *surf_fname)
		{
			printf("Writing surface to file %s...\n", surf_fname);

			// Open surface file
			FILE *outfile;
			if((outfile = fopen(surf_fname, "w")) == NULL) {
				printf("Error while trying to open %s for writing.\n", surf_fname);
				return -1;
			}

			fprintf(outfile, "surfacemesh\n");
			fprintf(outfile, "%d\n", num_nodes);
			for(int i = 0; i < num_nodes; i++) {
				fprintf(outfile, "%f %f %f\n", nodes[i].x, nodes[i].y, nodes[i].z);
			}

			fprintf(outfile, "%d\n", num_faces);

			for(std::list<Face*>::iterator face_iter = faces.begin() ; face_iter != faces.end(); face_iter++) {
				fprintf(outfile, "%d %d %d\n", (*face_iter)->get_node(1) + 1, (*face_iter)->get_node(2) + 1, (*face_iter)->get_node(3) + 1);
			}
			printf("...Done.\n");

			fclose(outfile);
			return 0;
		}

		void clean()
		{
			delete[] nodes;
			delete[] node_map;
			faces.clear();

			num_nodes = 0;
			num_nodes_original = 0;
			num_faces = 0;
			num_faces_original = 0;
			nodes = NULL;
			node_map = NULL;
		}

		private:
			int num_nodes_original, num_nodes, num_faces_original, num_faces;
			int *node_map;
			vector3 *nodes;
			std::list<Face*> faces, neighbourhood;
			
};

int main(int argc, char **argv)
{
	if(argc != 4) {
		printf("Usage: ./surface_coarse_grainer [SURF INPUT FNAME] [SURF OUTPUT FNAME] [LENGTH THRESHOLD]\n");
		return -1;
	}

	
	Surface surf;
	surf.init(argv[1]);
	surf.coarsen(atof(argv[3]));
	surf.output(argv[2]);
	surf.clean();
	return 0;
}
