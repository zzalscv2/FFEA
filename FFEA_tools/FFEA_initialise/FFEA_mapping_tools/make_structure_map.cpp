#include <cmath>
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include <string>
#include <string.h>
#include "Vectors.hpp"
#include "Matrices.hpp"
#define MAX_FNAME_SIZE 255

using namespace std;

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

		// Data members
		int num_nodes;
		int *n;
		vector3 **node;
};

int get_num_atoms_from_file(string pdb_fname) {

	// Open file
	FILE *fin;
	char buf[100];
	int num_atoms;
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

tet_element * find_containing_element(vector3 node, tet_element *top, int num_elements) {
	
	int i, j, num_faces_behind, check;
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
			return &top[i];
		}
	}

	// If we didn't trigger return condition, node is outside element
	return NULL;
}

tet_element * get_nearest_element(vector3 node, tet_element *top, int num_elements) {
	
	int i, j;
	double shortest_distance;
	vector3 distance;
	tet_element *nearest_element;
	shortest_distance = INFINITY;
	for(i = 0; i < num_elements; ++i) {
		distance = top[i].get_centroid() - node;
		if(distance.r < shortest_distance) {
			shortest_distance = distance.r;
			nearest_element = &top[i];
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
	if(argc != 5 and argc != 6) {
		cout << "Usage " << argv[0] << " [INPUT Base FFEA .node] [INPUT Base FFEA .top] [INPUT Target FFEA .node / Atomic .pdb] [OUTPUT .map fname] [scaling_factor (.pdb mapping only)]" << endl;
		exit(0);
	}

	// Get arguments
	string from_node_fname, to_fname, from_top_fname, output_map_fname;
	double scale;
	from_node_fname = argv[1];
	from_top_fname = argv[2];
	to_fname = argv[3];	
	output_map_fname = argv[4];
	if(argc == 6) {
		scale = atof(argv[5]);
	} else {
		scale = 1.0;
	}

	// Create node lists and topologies and map
	int num_elements_from, num_nodes_from, num_nodes_to;
	double **map;
	vector3 *from_node, *to_node;
	tet_element *from_top;
	
	// Sort out base structure
	cout << "Building base structure..." << endl << flush;
	num_nodes_from = get_num_nodes_from_file(from_node_fname);
	from_node = new vector3[num_nodes_from + 1];
	extract_and_create_nodes(from_node_fname, from_node);

	num_elements_from = get_num_elements_from_file(from_top_fname);
	from_top = new tet_element[num_elements_from];
	extract_and_create_topologies(from_top_fname, from_top, from_node);
	cout << "done!" << endl;

	// Now target structure, wich depends on the input type
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
		extract_and_create_pdb(to_fname, to_node, 1.0);
	}
	cout << "done!" << endl;

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
	int i, nearest_node_index, check;
	tet_element *containing_element, *temp_element;
	double single_map[4];
	vector3 nearest_node, translation;

	if(num_nodes_to > 20) {
		check = num_nodes_to / 20;
	} else {
		check = 1;
	}
	for(i = 0; i < num_nodes_to; ++i) {

		// Output check
		cout << "\r" << ((100 * i) / num_nodes_to) << "\% of map built...";

		// Find containing element
		containing_element = find_containing_element(to_node[i], from_top, num_elements_from);
		if(containing_element == NULL) {
			
			// Get nearest node from nearest element
			containing_element = get_nearest_element(to_node[i], from_top, num_elements_from);
		} 

		// Get position of to_node[i] as a function of the nodes of this element
		get_single_node_map(to_node[i], containing_element, single_map);

		// Add this data to final map
		add_map_to_big_map(single_map, map, i, containing_element);
		
	}
	cout << "done" << endl;

	// Print out the whole map
	cout << "Writing map to " << output_map_fname << "...";
	print_map_to_file(output_map_fname, map, num_nodes_from, num_nodes_to, 0);
	cout << "done! " << endl;
	return 0;
}
