#include <cmath>
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include <string>
#include <string.h>
#include <vector>
#include "Vectors.hpp"
#include "Matrices.hpp"
#define MAX_FNAME_SIZE 255

using namespace std;

struct distindex
{
	int index;
	double distance;
};

// Class for storing node connectivities (FFEA elements)
class tet_element
{
	public:
		
		// Constructors/Destructors
		tet_element() {
			num_nodes = 4;
			n = new int[num_nodes];
			node = new vector3*[num_nodes];
			for(int i = 0; i < num_nodes; ++i) {
				n[i] = 0;
				node[i] = NULL;
			}
		}
		
		tet_element(int num_nodes) {
			this->num_nodes = num_nodes;
			n = new int[num_nodes];
			node = new vector3*[num_nodes];
			for(int i = 0; i < num_nodes; ++i) {
				n[i] = 0;
				node[i] = NULL;
			}
		}
		
		tet_element(int n0, int n1, int n2, int n3, vector3 *nodes) {
			num_nodes = 4;
			n = new int[num_nodes];
			node = new vector3*[num_nodes];
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;
			n[3] = n3;
			node[0] = &nodes[n[0]];
			node[1] = &nodes[n[1]];
			node[2] = &nodes[n[2]];
			node[3] = &nodes[n[3]];
		}
		
		tet_element(int n0, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, vector3 *nodes) {
			num_nodes = 10;
			n = new int[num_nodes];
			node = new vector3*[num_nodes];
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;
			n[3] = n3;
			n[4] = n4;
			n[5] = n5;
			n[6] = n6;
			n[7] = n7;
			n[8] = n8;
			n[9] = n9;
			node[0] = &nodes[n[0]];
			node[1] = &nodes[n[1]];
			node[2] = &nodes[n[2]];
			node[3] = &nodes[n[3]];
			node[4] = &nodes[n[4]];
			node[5] = &nodes[n[5]];
			node[6] = &nodes[n[6]];
			node[7] = &nodes[n[7]];
			node[8] = &nodes[n[8]];
			node[9] = &nodes[n[9]];
		}
		
		~tet_element() {
			for(int i = 0; i < num_nodes; ++i) {
				n[i] = 0;
				node[i] = NULL;
			}
		}
		
		// Member functions
		void set_nodes_linear(int n0, int n1, int n2, int n3, vector3 *nodes) {
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;
			n[3] = n3;
			node[0] = &nodes[n[0]];
			node[1] = &nodes[n[1]];
			node[2] = &nodes[n[2]];
			node[3] = &nodes[n[3]];
		}
		
		vector3 get_face_normal(int index) {
			int n[3];			
			vector3 edge[2], to_centroid, normal;

			// Get face
			if(index == 0) {
				n[0] = 0;
				n[1] = 1;
				n[2] = 2;
				
			} else if(index == 1) {
				n[0] = 0;
				n[1] = 1;
				n[2] = 3;
				
			} else if(index == 2) {
				n[0] = 0;
				n[1] = 2;
				n[2] = 3;
				
			} else if(index == 3) {
				n[0] = 1;
				n[1] = 2;
				n[2] = 3;
				
			} else {
				cout << "Error. Invalid index, " << index << endl;
				exit(0);
			}

			// Get edge vectors
			edge[0] = *node[n[1]] - *node[n[0]];
			edge[1] = *node[n[2]] - *node[n[0]];

			// Get a normal
			normal = edge[0].cross(edge[1]);
			normal.normalise();

			// Make normal point out using projection onto centroid - n[0]
			to_centroid = get_centroid();
			to_centroid -= *node[n[0]];
			
			if(normal.dot(to_centroid) > 0.0) {
				normal *= -1;
			}
			return normal;
		}
		
		vector3 get_face_centroid(int index) {
			vector3 midpoint;
			if(index == 0) {
				midpoint = *node[0] + *node[1] + *node[2];
				midpoint *= 1.0 / 3.0;
			} else if(index == 1) {
				midpoint = *node[1] + *node[2] + *node[3];
				midpoint *= 1.0 / 3.0;
			} else if(index == 2) {
				midpoint = *node[0] + *node[2] + *node[3];
				midpoint *= 1.0 / 3.0;
			} else if(index == 3) {
				midpoint = *node[0] + *node[1] + *node[3];
				midpoint *= 1.0 / 3.0;
			} else {
				cout << "Error. Invalid index, " << index << endl;
				exit(0);
			}
			
			return midpoint;
		}
		
		vector3 get_centroid() {
			vector3 centroid = *node[0] + *node[1] + *node[2] + *node[3];
			centroid *= 0.25;
			return centroid;
		}

		int get_nearest_node_index(vector3 anode) {
			
			int nearest_node;
			double shortest_distance;
			vector3 distance;
			
			shortest_distance = INFINITY;
			for(int i = 0; i < 4; ++i) {
				distance = *node[i] - anode;
				if (distance.r < shortest_distance) {
					shortest_distance = distance.r;
					nearest_node = i;
				}
			}
			return nearest_node;
		}

		void print_detail() {
			cout << "Nodes: " << endl;
			for(int i = 0; i < 4;++i) {
				cout << "    " << n[i] << " (" << (*node[i]).x << ", " << (*node[i]).y << ", " << (*node[i]).z << ")" << endl;
			}
			cout << endl;
		}
		
		// Data members
		int num_nodes;
		int *n;
		vector3 **node;
		
};

int get_num_atoms_from_file(string pdb_fname) {

	// Open file
	FILE *fin;
	char buf[100];
	int num_atoms = 0;
	fin = fopen(pdb_fname.c_str(), "r");
	while(true) {
		if(feof(fin)) {
			break;	
		} else {
			fgets(buf,5,fin);
			if(strcmp(buf, "ATOM") == 0) {
				num_atoms++;
			}
			fgets(buf,96,fin);
		}
	}
	fclose(fin);
	cout << "\tFrom " + pdb_fname + " we expect " << num_atoms << " atoms." << endl;
	return num_atoms;
}

int get_num_nodes_from_file(string node_fname) {
	
	// Open file
	ifstream node_file;
	node_file.open(node_fname.c_str());
	
	string line;
	
	// Check correct file type
	getline(node_file, line);
	if (line.compare("ffea node file") != 0) {
		cout << "Error. Expected 'ffea node file', got '" << line << "'. Possibly not an ffea node file." << endl;
		exit(0);
	}
	
	// Get num_nodes
	int num_nodes;
	node_file >> line >> num_nodes;
	cout << "\tFrom " + node_fname + " we expect " << num_nodes << " nodes." << endl;
	return num_nodes;
}

int get_num_elements_from_file(string top_fname) {
	
	// Open file
	ifstream top_file;
	top_file.open(top_fname.c_str());
	
	string line;
	
	// Check correct file type
	getline(top_file, line);
	if (line.compare("ffea topology file") != 0) {
		cout << "Error. Expected 'ffea topology file', got '" << line << "'. Possibly not an ffea topology file." << endl;
		exit(0);
	}
	
	// Get num_elements
	int num_elements;
	top_file >> line >> num_elements;
	
	top_file.close();
	cout << "\tFrom " + top_fname + " we expect " << num_elements << " elements." << endl;
	return num_elements;
}

void extract_and_create_nodes(string node_fname, vector3 *node) {

	// Open file
	ifstream node_file;
	node_file.open(node_fname.c_str());
	
	string line;
	vector3 centroid;

	// Check correct file type
	getline(node_file, line);
	if (line.compare("ffea node file") != 0) {
		cout << "Error. Expected 'ffea node file', got '" << line << "'. Possibly not an ffea node file." << endl;
		exit(0);
	}
	
	// Get num_nodes and build node_list
	int num_nodes, num_surface_nodes, num_interior_nodes;
	node_file >> line >> num_nodes;
	cout << "\tFound " << num_nodes << " nodes." << endl;
	node_file >> line >> num_surface_nodes;
	node_file >> line >> num_interior_nodes;

	// Read in surface nodes
	int i;
	double x, y, z;
	getline(node_file, line);
	getline(node_file, line);
	for(i = 0; i < num_surface_nodes; ++i) {
		node_file >> x >> y >> z;
		node[i].set_pos(x, y, z);
		centroid.x += x;
		centroid.y += y;
		centroid.z += z;
	}
	
	// Read in interior nodes
	getline(node_file, line);
	getline(node_file, line);
	for(i = 0; i < num_interior_nodes; ++i) {
		node_file >> x >> y >> z;
		node[i + num_surface_nodes].set_pos(x, y, z);
		centroid.x += x;
		centroid.y += y;
		centroid.z += z;
	}
	
	// Close and return pointer
	centroid.x /= num_nodes;
	centroid.y /= num_nodes;
	centroid.z /= num_nodes;

	cout << "\tMoving to origin to overlap centroids...";
	for(i = 0; i < num_nodes; ++i) {
		node[i].x -= centroid.x;
		node[i].y -= centroid.y;
		node[i].z -= centroid.z;
	}
	cout << "done!" << endl;
	node_file.close();
	return;
}

void extract_and_create_pdb(string pdb_fname, vector3 *node, double scale) {

	// Open file
	FILE *fin;
	char buf[100];
	double x, y, z;
	fin = fopen(pdb_fname.c_str(), "r");
	int i = -1;
	vector3 centroid;
	while(true) {
		if(feof(fin)) {
			break;	
		} else {
			fgets(buf,5,fin);
			if(strcmp(buf, "END") == 0) {
				break;
			} else if(strcmp(buf, "ATOM") == 0) {
				i++;
				fgets(buf,27,fin);
				fgets(buf,9,fin);
				x = atof(buf) * scale;
				centroid.x += x;
				fgets(buf,9,fin);
				y = atof(buf) * scale;
				centroid.y += y;
				fgets(buf,9,fin);
				z = atof(buf) * scale;
				centroid.z += z;
				fgets(buf,100,fin);
				node[i].set_pos(x, y, z);
			} else {
				fgets(buf,100,fin);
			}
		}
	}
	int num_atoms = i + 1;
	cout << "\tFound " << num_atoms << " atoms." << endl;
	centroid.x /= num_atoms;
	centroid.y /= num_atoms;
	centroid.z /= num_atoms;

	cout << "\tMoving to origin to overlap centroids...";
	for(i = 0; i < num_atoms; ++i) {
		node[i].x -= centroid.x;
		node[i].y -= centroid.y;
		node[i].z -= centroid.z;
	}
	cout << "done!" << endl;
	fclose(fin);
	return;
}

void extract_and_create_topologies(string top_fname, tet_element *elem, vector3 *node) {

	// Open file
	ifstream top_file;
	top_file.open(top_fname.c_str());
	
	string line;
	
	// Check correct file type
	getline(top_file, line);
	if (line.compare("ffea topology file") != 0) {
		cout << "Error. Expected 'ffea topology file', got '" << line << "'. Possibly not an ffea topology file." << endl;
		exit(0);
	}
	
	// Get num_elements and build node_list
	int num_elements, num_surface_elements, num_interior_elements;
	top_file >> line >> num_elements;
	cout << "\tFound " << num_elements << " elements." << endl;
	top_file >> line >> num_surface_elements;
	top_file >> line >> num_interior_elements;

	// Read in surface elements
	int i;
	int n[10];
	getline(top_file, line);
	getline(top_file, line);
	for(i = 0; i < num_surface_elements; ++i) {
		top_file >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5] >> n[6] >> n[7] >> n[8] >> n[9];
		elem[i].set_nodes_linear(n[0], n[1], n[2], n[3], node);
	}
	
	// Read in interior elements
	getline(top_file, line);
	getline(top_file, line);
	for(i = 0; i < num_interior_elements; ++i) {
		top_file >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5] >> n[6] >> n[7] >> n[8] >> n[9];
		elem[i + num_surface_elements].set_nodes_linear(n[0], n[1], n[2], n[3], node);
	}
	
	// Close and return pointer
	top_file.close();
	return;
}

void build_pdb_topology(vector3 *from_node, vector3 *to_node, tet_element *from_top, int num_nodes_to, int num_nodes_from, int elements_per_node) {

	cout << "Building a pseudo-topology for the pdb structure..." << endl;
	
	// Now, to build these elements, we need 4 times that number of atoms
	int num_atoms_required = 4 * elements_per_node;
	
	// For every node requiring a map
	int i, j, k;
	distindex distance[num_atoms_required + 1], tempd;
	vector3 tempv, edge[3], normal;
	for(i = 0; i < num_nodes_to; ++i) {
	
		// Get a list of the closest 'num_atoms_required' atoms
		// This is lazy programming because I am busy! For faster code, please thing of a better way of doing this...
		cout << "\r" << ((100 * i) / num_nodes_to) << "\% structure built..." << flush;
		
		// Initialise array
		for(j = 0; j < num_atoms_required; ++j) {
			distance[j].index = -1;
			distance[j].distance = INFINITY;
		}
		
		// Get a new distance
		for(j = 0; j < num_nodes_from; ++j) {
			tempv = to_node[i] - from_node[j];
			distance[num_atoms_required].distance = tempv.r;
			//cout << to_node[i].r << " " << from_node[j].r << " " << tempv.r << endl;
			distance[num_atoms_required].index = j;
			
			// Insertion sort
			for(k = num_atoms_required; k > 0; k--) {
				if(distance[k].distance < distance[k - 1].distance) {
					tempd = distance[k];
					distance[k] = distance[k - 1];
					distance[k - 1] = tempd; 
				}
			}
		}
		//for(j = 0; j < num_atoms_required;++j) {
			//cout << distance[j].index << " " << distance[j].distance << endl;
		///	cout << "   " << to_node[i].x - from_node[distance[j].index].x << " " << to_node[i].y - from_node[distance[j].index].y << " " << to_node[i].z - from_node[distance[j].index].z << endl;
		//}
		//exit(0);
		
		// We need only the first 'num_atoms_required' indices to build elements out of!
		for(j = 0; j < elements_per_node; ++j) {
		
			// Now, build an element from this...
			// Get a face
			for(k = 0; k < 3; ++k) {
				edge[k] = from_node[distance[4 * j + k + 1].index] - from_node[distance[4 * j].index];
			}
		
			normal = edge[0].cross(edge[1]);
		
			// If this normal points positively projects into edge[3], then the normal points inwards. Reverse n[1] and n[2]
			//cout << distance[4 * j + k].index << " " << distance[4 * j + k + 1].index << " " << distance[4 * j + k + 2].index << " " << distance[4 * j + k + 3].index << endl;
			cout << endl;
			cout << to_node[i].x << " " << to_node[i].y << " " << to_node[i].z << endl;
			if(edge[2].dot(normal) > 0) {
				from_top[elements_per_node * i + j].set_nodes_linear(distance[0].index, distance[2].index, distance[1].index, distance[3].index, from_node);
				cout << from_node[distance[0].index].x << " " << from_node[distance[0].index].y << " " << from_node[distance[0].index].z << endl;
				cout << from_node[distance[2].index].x << " " << from_node[distance[2].index].y << " " << from_node[distance[2].index].z << endl; 
				cout << from_node[distance[1].index].x << " " << from_node[distance[1].index].y << " " << from_node[distance[1].index].z << endl; 
				cout << from_node[distance[3].index].x << " " << from_node[distance[3].index].y << " " << from_node[distance[3].index].z << endl; 
			} else {
				from_top[elements_per_node * i + j].set_nodes_linear(distance[0].index, distance[1].index, distance[2].index, distance[3].index, from_node);
				cout << from_node[distance[0].index].x << " " << from_node[distance[0].index].y << " " << from_node[distance[0].index].z << endl;
				cout << from_node[distance[1].index].x << " " << from_node[distance[1].index].y << " " << from_node[distance[1].index].z << endl; 
				cout << from_node[distance[2].index].x << " " << from_node[distance[2].index].y << " " << from_node[distance[2].index].z << endl; 
				cout << from_node[distance[3].index].x << " " << from_node[distance[3].index].y << " " << from_node[distance[3].index].z << endl;
			}
			exit(0);
		}
	}
	return;
}

/*void build_pdb_topology(vector3 *from_node, vector3 *to_node, tet_element *from_top, int num_nodes_to, int num_nodes_from) {
	
	cout << "Building a pseudo-topology for the pdb structure..." << endl;
	// For every node in to_node, find the closest 4 points in from_node and make that an element
	int i, j, k, n[4], tempn;
	vector3 d[4], tempd, distance, normal, edge[3];
	for (i = 0; i < num_nodes_to; ++i) {
		
		cout << "\r" << ((100 * i) / num_nodes_to) << "\% structure built...";
		
		// Initialise the element to the first 4 nodes
		for(j = 0; j < 4; ++j) {
			n[j] = j;
			d[j].r = INFINITY;
		}

		for(j = 3; j >=0; --j) {
			n[3] = j;
			d[3] = to_node[i] - from_node[j];

			for(k = 3; k >=1; --k) {
				if(d[k].r < d[k - 1].r) {
					tempd = d[k];
					d[k] = d[k - 1];
					d[k - 1] = tempd;
					tempn = n[k];
					n[k] = n[k - 1];
					n[k - 1] = tempn;
				} else {
					break;
				}
			}
		}
		
		for(j = 4; j < num_nodes_from; ++j) {
		
			// Get the distance to the target node from pdb node
			distance = to_node[i] - from_node[j];
			
			// Stick it in the element list
			if(distance.r < d[3].r) {
				d[3] = distance;
				n[3] = j;
			}
			
			for(k = 3; k >= 1; --k) {
				if(d[k].r < d[k - 1].r) {
					tempd = d[k];
					d[k] = d[k - 1];
					d[k - 1] = tempd;
					tempn = n[k];
					n[k] = n[k - 1];
					n[k - 1] = tempn;
				} else {
					break;
				}
			}
		}
		
		// Now, build an element from this...
		// Get a face
		for(j = 0; j < 3; ++j) {
			edge[j] = from_node[n[j + 1]] - from_node[n[0]];
		}
		
		normal = edge[0].cross(edge[1]);
		
		// If this normal points positively projects into edge[3], then the normal points inwards. Reverse n[1] and n[2]
		if(edge[2].dot(normal) > 0) {
			tempn = n[1];
			n[1] = n[2];
			n[2] = tempn;
		}
		
		from_top[i].set_nodes_linear(n[0], n[1], n[2], n[3], from_node);
	}
	cout << "\r" << ((100 * i) / num_nodes_to) << "\% structure built...";
}*/
		
void print_map_to_file(string map_fname, double **map, int num_nodes_from, int num_nodes_to, int sparsity) {

	int i, j;

	// Open file
	ofstream map_file;
	map_file.open(map_fname.c_str());

	// Write initial stuff
	map_file << "FFEA Kinetic Conformation Mapping File ";
	if(sparsity == 1) {
		map_file << "(Sparse)" << endl;
	} else {
		map_file << "(Dense)" << endl;
	}

	// Get num entries
	int num_entries = 0;
	for(i = 0; i < num_nodes_to; ++i) {
		for(j = 0; j < num_nodes_from; ++j) {

			// If irrelevant, don't include
			if (fabs(map[i][j]) >= 1e-5) {
				num_entries++;
			}
		}
	}
	map_file << "num_nodes_from " << num_nodes_from << endl;
	map_file << "num_nodes_to " << num_nodes_to << endl;
	map_file << "num_entries " << num_entries << endl;
	map_file << "map:" << endl;
	
	for(i = 0; i < num_nodes_to; ++i) {
		for(j = 0; j < num_nodes_from; ++j) {

			// If irrelevant, don't include
			if (fabs(map[i][j]) < 1e-5) {
				map_file << 0 << " ";
			} else {
				map_file << map[i][j] << " ";
			}
		}
		map_file << endl;
	}
	return;
}

tet_element * find_containing_element(vector3 node, tet_element *top, int num_elements, vector<tet_element*> used) {
	
	bool failed;
	int i, j, num_faces_behind, check;
	vector<tet_element*>::iterator it;
	//double angle;
	vector3 face_normal, face_centroid, face_to_node;
	check = 0;

	for(i = 0; i < num_elements; ++i) {

		// Check if in element by calculating face normals and node-face_center vectors
		//cout << "Element 0: " << top[i].n[0] << " " << top[i].n[1] << " " << top[i].n[2] << " " << top[i].n[3] << endl;
		num_faces_behind = 0;
		for(j = 0; j < 4; ++j) {

			face_normal = top[i].get_face_normal(j);
			//cout << "\tFace normal " << j << ": " << face_normal.x << " " << face_normal.y << " " << face_normal.z << endl;
			face_centroid = top[i].get_face_centroid(j);
			face_to_node = node - face_centroid;
			//cout << "\tFace to node " << j << ": " << face_to_node.x << " " << face_to_node.y << " " << face_to_node.z << endl;
			face_to_node.normalise();

			// If in front of plane, not in element. If behind all, it is in element
			if(face_normal.dot(face_to_node) > 0) {
				
				// Not in this element
				break;
			} else {
				num_faces_behind++;
			}
		}

		// If behind all faces, node is in element
		if(num_faces_behind == 4) {
		
			// Must check whether or not we have already used it
			failed = false;
			for(it = used.begin(); it != used.end(); ++it) {
				if(*it == &top[i]) {
					failed = true;
					break;
				}
			}
			if(!failed) {
				return &top[i];
			}
		}
	}

	// If we didn't trigger return condition, node is outside element
	return NULL;
}

tet_element * get_nearest_element(vector3 node, tet_element *top, int num_elements, vector<tet_element*> used) {
	
	bool failed;
	int i, j;
	double shortest_distance;
	vector3 distance;
	tet_element *nearest_element;
	vector<tet_element*>::iterator it;
	shortest_distance = INFINITY;
	for(i = 0; i < num_elements; ++i) {
		distance = top[i].get_centroid() - node;
		if(distance.r < shortest_distance) {
		
			// Check we haven't used it
			failed = false;
			for(it = used.begin(); it != used.end(); ++it) {
				if(*it == &top[i]) {
					failed = true;
					break;
				}
			}
			if(!failed) {
				shortest_distance = distance.r;
				nearest_element = &top[i];
			}
		}
	}

	return nearest_element;
}

void get_single_node_map(vector3 node, tet_element *elem, double *map) {

	vector3 a, b, c, n, coeff;
	matrix33 M, M_inv;

	// Get the required vectors and make them a matrix
	a = *elem->node[1] - *elem->node[0];
	b = *elem->node[2] - *elem->node[0];
	c = *elem->node[3] - *elem->node[0];
	n = node - *elem->node[0];
	M.set(a, b, c);
	M_inv = M.get_inverse();
	coeff = M_inv.apply(n);
	map[0] = 1 - coeff.x - coeff.y - coeff.z;
	map[1] = coeff.x;
	map[2] = coeff.y;
	map[3] = coeff.z;

	return;
}

tet_element * make_containing_element(tet_element *old_containing_element, int nearest_node_index, int new_node_index, vector3 *from_node) {

	int i, new_node[4];
	tet_element *new_containing_element;
	new_containing_element = new tet_element;

	// Get and set new nodes
	for(i = 0; i < 4; ++i) {
		if(i == nearest_node_index) {
			new_node[i] = new_node_index;
		} else {
			new_node[i] = old_containing_element->n[i];
		}
	}
	new_containing_element->set_nodes_linear(new_node[0], new_node[1], new_node[2], new_node[3], from_node);
	return new_containing_element;
}

void add_map_to_big_map(double *little_map, double **big_map, int to, tet_element *containing_elem) {

	int i;
	for(i = 0; i < 4; ++i) {
		big_map[to][containing_elem->n[i]] = little_map[i];
	}
	return;
}

int main(int argc, char **argv) {

	// Check for the right input
	if((argc - 1 )% 2 != 0 || argc < 4) {
		cout << "Usage " << argv[0] << " -i [INPUT Base FFEA .node/.pdb] -t [INPUT Base FFEA .top] -o [INPUT Target FFEA .node / Atomic .pdb] -m [OUTPUT .map fname] -s [pdb scaling_factor (.pdb mapping only)]" << endl;
		exit(0);
	}
	
	// Get arguments
	string from_node_fname, to_fname, from_top_fname, output_map_fname;
	double scale = 1.0;
	
	int set_in_pos = 0, set_in_top = 0, set_out_pos = 0, set_out_map = 0;
	
	// Get arguments in pairs
	string f1, f2, f3, f4, f5;
	f1 = "-i";
	f2 = "-t";
	f3 = "-o";
	f4 = "-m";
	f5 = "-s"; 
	for(int i = 1; i < argc; i += 2) {
		if (f1.compare(argv[i]) == 0) {
			from_node_fname = argv[i + 1];
			set_in_pos = 1;
		} else if (f2.compare(argv[i]) == 0) {
			from_top_fname = argv[i + 1];
			set_in_top = 1;
		}  else if (f3.compare(argv[i]) == 0) {
			to_fname = argv[i + 1];
			set_out_pos = 1;
		}  else if (f4.compare(argv[i]) == 0) {
			output_map_fname = argv[i + 1];
			set_out_map = 1;
		}  else if (f5.compare(argv[i]) == 0) {
			scale = atof(argv[i + 1]);
		} else {
			cout << "Unrecognised flag" << endl;
			exit(0);
		}
	}

	// Create node lists and topologies and map
	int num_elements_from, num_nodes_from, num_nodes_to;
	int elements_per_node;
	double **map;
	vector3 *from_node, *to_node;
	tet_element *from_top;
	
	// Build target structure first, which depends on the input type
	cout << "Building target structure..." << endl << flush;
	string line;
	ifstream fin;
	fin.open(to_fname.c_str());
	getline(fin, line);
	fin.close();

	if(line == "ffea node file" || line == "ffea node file\n") {
		num_nodes_to = get_num_nodes_from_file(to_fname);
		to_node = new vector3[num_nodes_to];
		extract_and_create_nodes(to_fname, to_node);
	} else {
		num_nodes_to = get_num_atoms_from_file(to_fname);
		to_node = new vector3[num_nodes_to];
		extract_and_create_pdb(to_fname, to_node, scale);
	}
	cout << "done!" << endl;
	
	// Now sort out base structure, which depends on input type
	cout << "Building base structure..." << endl << flush;
	fin.open(from_node_fname.c_str());
	getline(fin, line);
	fin.close();
	if(line == "ffea node file" || line == "ffea node file\n") {
	
		// No node target here!
		cout << "Error. Expected a pdb file" << endl;
		exit(0);
		
	} else {
	
		// Atomic positions
		num_nodes_from = get_num_atoms_from_file(from_node_fname);
		from_node = new vector3[num_nodes_from + 1];
		extract_and_create_pdb(from_node_fname, from_node, scale);
		
		// Build a topology structure using nearest nodes
			
		// There will be more atoms than nodes, by far. So, one element per node not enough!
		// Each element contains one unique atom (solid state was useful after all!) so atoms_per_node = num_elements_required
		elements_per_node = (num_nodes_from) / num_nodes_to;
		num_elements_from = num_nodes_to * elements_per_node;
		
		from_top = new tet_element[num_elements_from];
		build_pdb_topology(from_node, to_node, from_top, num_nodes_to, num_nodes_from, elements_per_node);			// Fix me		
	}
	
	

	// Build the map object
	map = new double*[num_nodes_to];
	cout << "Building map object..." << endl << flush;
	for(int i = 0; i < num_nodes_to; ++i) {
		cout << "\r" << ((100 * i) / num_nodes_to) << "\% of map initialised...";
		map[i] = new double[num_nodes_from];
		for(int j = 0; j < num_nodes_from; ++j) {
			map[i][j] = 0.0;
		}
	}
	cout << "done!" << endl;
	
	// Find nearest LINEAR node, and therefore nearest element / containing element
	int i, j, nearest_node_index;
	tet_element *containing_element, *temp_element;
	double single_map[4];
	vector3 nearest_node, translation, partial_node;
	vector<tet_element*> used;
	for(i = 0; i < num_nodes_to; ++i) {

		// Reset list of used elements
		used.clear();
		
		// Output check
		cout << "\r" << ((100.0 * i) / num_nodes_to) << "\% of map built..." << flush;

		// Now, we must use 'elements_per_node' elements
		
		// Each element contributes equally to the position of the new node
		partial_node = to_node[i] * (1.0 / elements_per_node);
		for(j = 0; j < elements_per_node; ++j) {
		
			// Find containing element
			containing_element = find_containing_element(to_node[i], from_top, num_elements_from, used);
			if(containing_element == NULL) {
			
				// Get nearest node from nearest element
				containing_element = get_nearest_element(to_node[i], from_top, num_elements_from, used);
			} 

			//containing_element->print_detail();
			used.push_back(containing_element);
			
			// Get position of to_node[i] as a function of the nodes of this element
			//get_single_node_map(to_node[i], containing_element, single_map);
			get_single_node_map(partial_node, containing_element, single_map);
			
			// Add this data to final map
			add_map_to_big_map(single_map, map, i, containing_element);
		}
	}
	cout << "done" << endl;

	// Print out the whole map
	cout << "Writing map to " << output_map_fname << "...";
	print_map_to_file(output_map_fname, map, num_nodes_from, num_nodes_to, 0);
	cout << "done! " << endl;
	return 0;
}
