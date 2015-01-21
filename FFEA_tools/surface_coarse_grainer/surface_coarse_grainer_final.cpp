#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <string.h>
#include <unistd.h>

using namespace std;
class vector3
{
	public:
		
		vector3() 
		{
			this->x = 0.0;
			this->y = 0.0;
			this->z = 0.0;
		}
		
		vector3(double x, double y, double z) 
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}
		
		// Vector functions
		void normalise() {
			double l = sqrt(x * x + y * y + z * z);
			x /= l;
			y /= l;
			z /= l;
		}
		
		vector3 get_normal() {
			double l = sqrt(x * x + y * y + z * z);
			return vector3(this->x/l, this->y/l, this->z/l);
		}
		
		double get_magnitude() {
			return sqrt(x * x + y * y + z * z);
		}
		
		void set_midpoint(vector3 a, vector3 b) {
			 this->x = (a.x + b.x) / 2.0;
			 this->y = (a.y + b.y) / 2.0;
			 this->z = (a.z + b.z) / 2.0;
		}
		
		void set_pos(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
			
		}
		
		void cross_prod(vector3 a, vector3 b) {
			this->x = a.y * b.z + b.y * a.z;
			this->y = -1 * (a.x * b.z - a.z * b.x);
			this->z = a.x * b.y - a.y * b.z; 	
		}
		
		// Vector components
		double x;
		double y;
		double z;
};

class Edge
{
	public:
	
		Edge()
		{
			this->n_index[0] = 0;
			this->n_index[1] = 0;
			this->length = 0.0;	
		}
		
		Edge(int n0, int n1, vector3 *node)
		{
			this->n_index[0] = n0;
			this->n_index[1] = n1;
			this->length = calc_length(node);
		}
		
		// Edge functions
		int init(int n0, int n1, vector3 *node) {
			this->n_index[0] = n0;
			this->n_index[1] = n1;
			this->length = calc_length(node);
		}
		double calc_length(vector3 *node) {
			vector3 v(node[this->n_index[1]].x - node[this->n_index[0]].x, node[this->n_index[1]].y - node[this->n_index[0]].y, node[this->n_index[1]].z - node[this->n_index[0]].z);
			this->length = v.get_magnitude();
			return v.get_magnitude();
		}
			
		// Edge components
		int n_index[2];
		double length;
};
class Face
{
	public:
		
		Face(int n0, int n1, int n2, vector3 *node)
		{
			this->n[0] = n0;
			this->n[1] = n1;
			this->n[2] = n2;
			this->remaining_node = n0;
			edge = new Edge[3];
			edge[0].init(n0, n1, node);
			edge[1].init(n1, n2, node);
			edge[2].init(n2, n0, node); 
			calc_shortest_edge();
		}
		
		// Face functions
		void calc_shortest_edge() {
			shortest_edge = INFINITY;
			for(int i = 0; i < 3; ++i) {
				if(edge[i].length < shortest_edge) {
					shortest_edge = edge[i].length;
					shortest_edge_index = i;
				}
			}
		}
		
		double get_edge_length(int n_a, int n_b) {
			for(int i = 0; i < 3; ++i) {
				if((edge[i].n_index[0] == n_a && edge[i].n_index[1] == n_b) || (edge[i].n_index[0] == n_b && edge[i].n_index[1] == n_a)) {
					return edge[i].length;
				}
			}
		}
		
		// Face components
		int n[3];
		int remaining_node;
		Edge *edge;
		double shortest_edge;
		int shortest_edge_index;
};

// Main class defining surface
class Surface
{
	public:
		
		Surface() 
		{
			node = NULL;
			num_nodes = 0;
			num_faces = 0;
			num_nodes_initial = 0;
			num_faces_initial = 0;
		}
		
		~Surface()
		{
			delete[] node;
			node = NULL;
			num_nodes = 0;
			num_faces = 0;
			num_nodes_initial = 0;
			num_faces_initial = 0;
		}
		
		// Surface functions
		int init(char *surf_fname) {
		
			int i, j;
			
			// Open file
			FILE *surf_file;
			surf_file = fopen(surf_fname, "r");
			if(surf_file == NULL) {
				printf("Error. Filename %s does not exist.\n", surf_fname);
				return -1;
			}
			
			// Check file type is correct
			char line[255];
			if(fscanf(surf_file, "%s\n", &line) != 1) {
				printf("Error. 'surfacemesh' line does not exist.\n");
				return -1;
			}
			if(strcmp(line, "surfacemesh") != 0) {
				printf("Error. File %s is not a NETGEN surface file\n");
				return -1;
			}
			
			// Get num_nodes
			if(fscanf(surf_file, "%d\n", &num_nodes) != 1) {
				printf("Error. 'num_nodes' line does not exist.\n");
				return -1;
			};
			num_nodes_initial = num_nodes;
			
			// Create nodes
			node = new vector3[num_nodes];
			
			// Read nodes
			for(i = 0; i < num_nodes; ++i) {
				if(fscanf(surf_file, "%lf %lf %lf\n", &node[i].x, &node[i].y, &node[i].z) != 3) {
					printf("Error. Expected %d nodes, found only %d.\n", num_nodes, i);
				}
			}
			
			// Get num_faces
			if(fscanf(surf_file, "%d\n", &num_faces) != 1) {
				printf("Error. 'num_faces' line does not exist.\n");
				return -1;
			};
			num_faces_initial = num_faces;
			
			// Read faces
			int n0, n1, n2;
			for(i = 0; i < num_faces; ++i) {
				if(fscanf(surf_file, "%d %d %d\n", &n0, &n1, &n2) != 3) {
					printf("Error. Expected %d faces, found only %d.\n", num_faces, i);
				}
				
				// Add to list
				face.push_back(new Face(n0 - 1, n1 - 1, n2 - 1, node));
			}
			
			// Close file and return 
			fclose(surf_file);
			return 0;
		}
		
		int coarsen(double limit, char *conserve_vol, char *find_smallest) {
			
			// Iterate through faces while any edge is less than limit
			int i, j, k, runs = 0, completed_check, num_shared_faces, num_shared_nodes;
			int node_to_delete[2], remaining_node, do_not_delete, section_deleted, twod_section_exists, face_sharing_remaining_node[2];
			double length_deleted;
			vector3 midpoint;
			list<Face*>::iterator face_iterator, face_to_delete[4];
			while(true) {
				
				// Output completion data
				// Output Info
				if((num_faces_initial - num_faces) % 1 == 0) {
					printf("Faces Remaining = %d\n", num_faces);
					printf("Nodes Remaining = %d\n", num_nodes);
					printf("Last Length Deleted = %f\n\n", length_deleted);					
				}

				runs++;
				completed_check = 1;

				// Find an edge requiring deletion, flag face and two nodes
				if(strcmp(find_smallest, "y") == 0 || strcmp(find_smallest, "Y") == 0) {

					// Find smallest edge in system first and face containing it
					length_deleted = INFINITY;
					for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
						if((*face_iterator)->shortest_edge < length_deleted) {
							length_deleted = (*face_iterator)->shortest_edge;
							face_to_delete[0] = face_iterator;
						}
					}

					// Not finished if length < limit
					if(length_deleted < limit) {
						completed_check = 0;
					}
		
					// Get nodes to delete
					for(i = 0; i < 2; ++i) {
						node_to_delete[i] = (*face_to_delete[0])->edge[(*face_to_delete[0])->shortest_edge_index].n_index[i];
					}

					// Get remaining node in that face
					for(i = 0; i < 3; ++i) {
						if((*face_to_delete[0])->n[i] != node_to_delete[0] && (*face_to_delete[0])->n[i] != node_to_delete[1]) {
							(*face_to_delete[0])->remaining_node = (*face_to_delete[0])->n[i];
						}
					}
				} else {

					// Find any edge less than limit and face contatining it
					for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
						if((*face_iterator)->shortest_edge < limit) {
							length_deleted = (*face_iterator)->shortest_edge;
							face_to_delete[0] = face_iterator;
							
							// Get nodes to delete
							for(i = 0; i < 2; ++i) {
								node_to_delete[i] = (*face_iterator)->edge[(*face_iterator)->shortest_edge_index].n_index[i];
							}
							for(i = 0; i < 3; ++i) {
								if((*face_iterator)->n[i] != node_to_delete[0] && (*face_iterator)->n[i] != node_to_delete[1]) {
									(*face_iterator)->remaining_node = (*face_iterator)->n[i];
								}
							}
							completed_check = 0;
							break;
						}
					}
				}

				// Break loop if completed for remapping nodes to create a consistent surface
				if(completed_check == 1) {
					printf("\nCoarsening Completed.\n");
					break;
				}
				// Find the other face and the neighbourhood about the nodes
				// Empty previous neighbourhood
				neighbourhood.clear();

				// Find faces and neighbourhood based on how many of the nodes_to_delete they are connected to
				num_shared_faces = 0;
				for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
					
					num_shared_nodes = 0;
					for(i = 0; i < 3; ++i) {
						for(j = 0; j < 2; ++j) {
							if((*face_iterator)->n[i] == node_to_delete[j]) {
								num_shared_nodes++;
								break;
							}
						}
					}

					// If only shares one node, counts as a neighbourhood face
					if(num_shared_nodes == 1) {
						neighbourhood.push_back(*face_iterator);
					}

					// If shares two nodes, shares the whole edge. These faces need dealing with (may be more than two)
					// Will find original face, then overwrite with the second
					if(num_shared_nodes == 2) {
					
						// Remember how many faces are shared (2 is easy, 4 is not)
						num_shared_faces++;
			
						// Remember remaining node
						for(i = 0; i < 3; ++i) {
							if((*face_iterator)->n[i] != node_to_delete[0] && (*face_iterator)->n[i] != node_to_delete[1]) {
								(*face_iterator)->remaining_node = (*face_iterator)->n[i];
							}
						}
						face_to_delete[num_shared_faces - 1] = face_iterator;
					}
				}

				// Coarsen in different ways depending on how many faces are sharing the edge
				section_deleted = 0;
				//printf("Num_shared_faces = %d\n", num_shared_faces);
				if(num_shared_faces == 2) {
					// Move lowest node to midpoint
					midpoint.set_midpoint(node[node_to_delete[0]], node[node_to_delete[1]]);
					node[node_to_delete[0]].set_pos(midpoint.x, midpoint.y, midpoint.z);

					// Move highest node to infinity
					node[node_to_delete[1]].set_pos(INFINITY, INFINITY, INFINITY);
					num_nodes--;

					// Reconnect all neighbourhood faces to other node and recalculate shortest edge
					for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
					
						// Node indices
						for(i = 0; i < 3 ; ++i) {
							if((*face_iterator)->n[i] == node_to_delete[1]) {
								(*face_iterator)->n[i] = node_to_delete[0];
							}
						}	
					
						// Edge lengths
						for(i = 0; i < 3; ++i) {
							for(j = 0; j < 2; ++j) {
								if((*face_iterator)->edge[i].n_index[j] == node_to_delete[1]) {
									(*face_iterator)->edge[i].n_index[j] = node_to_delete[0];
									break;	
								}
							}
							(*face_iterator)->edge[i].calc_length(node);
						}
					
						// Shortest length
						(*face_iterator)->calc_shortest_edge();
					}

					// Delete faces
					face.erase(face_to_delete[1]);
					face.erase(face_to_delete[0]);
					num_faces -= 2;

				} else if(num_shared_faces == 4) {
					// Firstly check for a tetrahedral collapse (leaves two triangles sticking out of an otherwise consistent surface)
					// Delete the two evil triangles
				//	for(i = 0; i < num_shared_faces; ++i) {
				//		printf("Face %d: %d %d %d\n", i, (*face_to_delete[i])->n[0], (*face_to_delete[i])->n[1], (*face_to_delete[i])->n[2]);
				//	}

					twod_section_exists = 0;
					for(i = 0; i < num_shared_faces; ++i) {
						for(j = i + 1; j < num_shared_faces; ++j) {
							if((*face_to_delete[i])->remaining_node == (*face_to_delete[j])->remaining_node) {
								twod_section_exists = 1;
							//	printf("Face %d : %d %d %d\n", i, (*face_to_delete[i])->n[0], (*face_to_delete[i])->n[1], (*face_to_delete[i])->n[2]);
							//	printf("Face %d : %d %d %d\n", j, (*face_to_delete[j])->n[0], (*face_to_delete[j])->n[1], (*face_to_delete[j])->n[2]);
								face_sharing_remaining_node[0] = i;
								face_sharing_remaining_node[1] = j;
								remaining_node = (*face_to_delete[i])->remaining_node;
								
								// Exit these loops
								j = num_shared_faces;
								i = num_shared_faces;
								break;
							}
						}
					}
					
					if(twod_section_exists == 1) {
								
						// Erase the triangles
						face.erase(face_to_delete[face_sharing_remaining_node[0]]);
						face.erase(face_to_delete[face_sharing_remaining_node[1]]);
						num_faces -= 2;

						// Check that this node isn't contained anywhere else
						// If it is, we have created two separate surfaces, but they're at least consistent!
						do_not_delete = 0;
						for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
							for(k = 0; k < 3; ++k) {
								if((*face_iterator)->n[k] == remaining_node) {
									do_not_delete = 1;
					
									// Exit both loops
									face_iterator = face.end();
									face_iterator--;
									break;
								}				
							}
						}

						// If it isn't, delete it
						if(do_not_delete == 0) {
							printf("%d", remaining_node);
							node[remaining_node].set_pos(INFINITY, INFINITY, INFINITY);
							num_nodes--;
						}
						
					} else {

						// If not, then we have a very complicated structure
						// 4-way coarsening may work. If there is a seg fault, or a bug, it is likely here
						// Volume conservation may be possible via more accurate definition of the surface structure, 
						// but currently not implemented
						
						// Move lowest node to midpoint
						midpoint.set_midpoint(node[node_to_delete[0]], node[node_to_delete[1]]);
						node[node_to_delete[0]].set_pos(midpoint.x, midpoint.y, midpoint.z);
	
						// Move highest node to infinity
						node[node_to_delete[1]].set_pos(INFINITY, INFINITY, INFINITY);
						num_nodes--;
	
						// Reconnect all neighbourhood faces to other node and recalculate shortest edge
						for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
						
							// Node indices
							for(i = 0; i < 3 ; ++i) {
								if((*face_iterator)->n[i] == node_to_delete[1]) {
									(*face_iterator)->n[i] = node_to_delete[0];
								}
							}	
						
							// Edge lengths
							for(i = 0; i < 3; ++i) {
								for(j = 0; j < 2; ++j) {
									if((*face_iterator)->edge[i].n_index[j] == node_to_delete[1]) {
										(*face_iterator)->edge[i].n_index[j] = node_to_delete[0];
										break;	
									}
								}
								(*face_iterator)->edge[i].calc_length(node);
							}
						
							// Shortest length
							(*face_iterator)->calc_shortest_edge();
						}
	
						// Delete faces
						face.erase(face_to_delete[3]);
						face.erase(face_to_delete[2]);
						face.erase(face_to_delete[1]);
						face.erase(face_to_delete[0]);
						num_faces = num_faces - 4;
				
						printf("4-way coarsening attempted. Good luck!\n");
					}
				}
				
				// Check neighbourhood consistency
				for(face_iterator = neighbourhood.begin(); face_iterator != neighbourhood.end(); ++face_iterator) {
					
					// Checking that faces haven't collapsed
					for(i = 0; i < 3; ++i) {
						for(j = i + 1; j < 3; ++j) {
							if((*face_iterator)->n[i] == (*face_iterator)->n[j]) {
								printf("In neighbourhood checker\n");
								exit(0);
								face.erase(face_iterator);
								num_faces--;
							}
						}
					}
				}
			}
			
			// Remap all nodes
			node_map = new int[num_nodes_initial];
			int new_node_index = 0;
			for(i = 0; i < num_nodes_initial; ++i) {
				if(node[i].x != INFINITY) {
					node_map[i] = new_node_index;
					if(i == 433) {
						printf("Node %d -> %d\n", i, new_node_index + 1);
					}
					if(i == 278) {
						printf("Node %d -> %d\n", i, new_node_index + 1);
					}
					if(i == 280) {
						printf("Node %d -> %d\n", i, new_node_index + 1);
					}
					node[new_node_index] = node[i];
					new_node_index++;
				} else {
					node_map[i] = -1;
				}
			}
			
			// Remap nodes in faces
			for(face_iterator = face.begin(); face_iterator != face.end(); face_iterator++) {
				for(i = 0; i < 3; ++i) {
					(*face_iterator)->n[i] = node_map[(*face_iterator)->n[i]];
				}
			}
			return 0;
		}
		
		int write_to_file(char *surf_fname) {
			int i;
			FILE *surf_out;
			list<Face*>::iterator face_iterator;
			surf_out = fopen(surf_fname, "w");
			fprintf(surf_out, "surfacemesh\n");
			fprintf(surf_out, "%d\n", num_nodes);
			for(i = 0; i < num_nodes; ++i) {
				fprintf(surf_out, "%6.3lf %6.3lf %6.3lf\n", node[i].x, node[i].y, node[i].z);
			}
			fprintf(surf_out, "%d\n", num_faces);
			for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator) {
				fprintf(surf_out, "%d %d %d\n", (*face_iterator)->n[0] + 1, (*face_iterator)->n[1] + 1, (*face_iterator)->n[2] + 1);
			}
			fclose(surf_out);
			return 0;
		}
		
		void check_edges() {
			int i, j, index0, index1;
			list<Face*>::iterator face_iterator, face_iterator2;
			for(face_iterator = face.begin(); face_iterator != face.end(); ++face_iterator){
				for(i = 0; i < 3; ++i) {
					index0 = (*face_iterator)->edge[i].n_index[0];
					index1 = (*face_iterator)->edge[i].n_index[1];
					face_iterator2 = face_iterator;
					face_iterator2++;
					while(face_iterator2 != face.end()) {
						for(j = 0; j < 3; ++j) {
							if((*face_iterator2)->edge[j].n_index[0] == index0 && (*face_iterator2)->edge[j].n_index[1] == index1) {
								printf("Edge %d %d has length %lf\n", (*face_iterator)->edge[i].n_index[0], (*face_iterator)->edge[i].n_index[1], (*face_iterator)->edge[i].length);
								printf("Edge %d %d has length %lf\n\n", (*face_iterator2)->edge[j].n_index[0], (*face_iterator2)->edge[j].n_index[1], (*face_iterator2)->edge[j].length);
							} else if ((*face_iterator2)->edge[j].n_index[0] == index1 && (*face_iterator2)->edge[j].n_index[1] == index0) {
								printf("Edge %d %d has length %lf\n", (*face_iterator)->edge[i].n_index[0], (*face_iterator)->edge[i].n_index[1], (*face_iterator)->edge[i].length);
								printf("Edge %d %d has length %lf\n\n", (*face_iterator2)->edge[j].n_index[0], (*face_iterator2)->edge[j].n_index[1], (*face_iterator2)->edge[j].length);
							}
						}
						face_iterator2++;
					}
				}
			}
		}
	private:

		// Surface functions
		double calc_vol(list<Face*> local_surf, vector3 *node, vector3 normal, vector3 extra_node) {
			list<Face*>::iterator face_iterator;
			int i;
			
			exit(0);
			
		}
		
		// Surface components
		int num_nodes;
		int num_faces;
		int num_nodes_initial;
		int num_faces_initial;
		int *node_map;
		vector3 *node; 
		list<Face*> face, neighbourhood;
};

// Main program begins
int main(int argc, char **argv) {
	if(argc != 6) {
		cout << "Input Error" << endl << "USAGE: " << argv[0] << " [Input .surf fname] [Output .surf fname] [Coarseness level] [Volume conserve (y/n)] [Find smallest edge? (y/n)]" << endl;
		return -1;
	}
	
	// Create a surface
	Surface *surf = new Surface();
	if(surf->init(argv[1]) != 0) {
		printf("Error initialising surface\n.");
		return -1;
	};
	
	// Check Identical Edges are the same length
	//surf->check_edges();
	
	// Coarsen surface
	surf->coarsen(atof(argv[3]), argv[4], argv[5]);
	
	// Get rid of surface errors
	//surf->clean();
	
	// Output surface filebuf
	surf->write_to_file(argv[2]);
	
	// Delete surface
	delete surf;
	
	return 0;
}
