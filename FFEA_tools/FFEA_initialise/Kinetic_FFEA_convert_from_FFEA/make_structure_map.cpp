#include <cmath>
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include <string>
#include "Vectors.hpp"
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
			for(int i = 0; i < num_nodes; ++i) {
				n[i] = 0;
			}
		}
		
		tet_element(int num_nodes) {
			this->num_nodes = num_nodes;
			n = new int[num_nodes];
			for(int i = 0; i < num_nodes; ++i) {
				n[i] = 0;
			}
		}
		
		tet_element(int n0, int n1, int n2, int n3) {
			num_nodes = 4;
			n = new int[num_nodes];
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;
			n[3] = n3;
		}
		
		tet_element(int n0, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9) {
			num_nodes = 10;
			n = new int[num_nodes];
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
		}
		
		~tet_element() {
			for(int i = 0; i < num_nodes; ++i) {
				n[i] = 0;
			}
		}
		
		// Member functions
		void set_nodes_linear(int n0, int n1, int n2, int n3) {
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;
			n[3] = n3;
		}
		
		vector3 get_face_normal(int index, vector3 *nodes) {
			vector3 a, b, normal;
			if(index == 0) {
				a = nodes[n[1]] - nodes[n[0]];
				b = nodes[n[2]] - nodes[n[1]];
			} else if(index == 1) {
				a = nodes[n[2]] - nodes[n[3]];
				b = nodes[n[1]] - nodes[n[2]];
			} else if(index == 2) {
				a = nodes[n[3]] - nodes[n[2]];
				b = nodes[n[0]] - nodes[n[3]];
			} else if(index == 3) {
				a = nodes[n[0]] - nodes[n[1]];
				b = nodes[n[3]] - nodes[n[0]];
			} else {
				cout << "Error. Invalid index, " << index << endl;
				exit(0);
			}
			
			normal = a.cross(b);
			normal.normalise();
			return normal;
		}
		
		vector3 get_face_midpoint(int index, vector3 *nodes) {
			vector3 midpoint;
			if(index == 0) {
				midpoint = nodes[n[0]] + nodes[n[1]] + nodes[n[2]];
				midpoint *= 1.0 / 3.0;
			} else if(index == 1) {
				midpoint = nodes[n[1]] + nodes[n[2]] + nodes[n[3]];
				midpoint *= 1.0 / 3.0;
			} else if(index == 2) {
				midpoint = nodes[n[0]] + nodes[n[2]] + nodes[n[3]];
				midpoint *= 1.0 / 3.0;
			} else if(index == 3) {
				midpoint = nodes[n[0]] + nodes[n[1]] + nodes[n[3]];
				midpoint *= 1.0 / 3.0;
			} else {
				cout << "Error. Invalid index, " << index << endl;
				exit(0);
			}
			
			return midpoint;
		}
		
		// Data members
		int num_nodes;
		int *n;
};

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
	
	node_file.close();
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
	return num_elements;
}

void extract_and_create_nodes(string node_fname, vector3 *node) {

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
	
	// Get num_nodes and build node_list
	int num_nodes, num_surface_nodes, num_interior_nodes;
	node_file >> line >> num_nodes;
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
	}
	
	// Read in interior nodes
	getline(node_file, line);
	getline(node_file, line);
	for(i = 0; i < num_interior_nodes; ++i) {
		node_file >> x >> y >> z;
		node[i + num_surface_nodes].set_pos(x, y, z);
	}
	
	// Close and return pointer
	node_file.close();
	return;
}

void extract_and_create_topologies(string top_fname, tet_element *elem) {

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
	top_file >> line >> num_surface_elements;
	top_file >> line >> num_interior_elements;

	// Read in surface elements
	int i;
	int n0, n1, n2, n3;
	getline(top_file, line);
	getline(top_file, line);
	for(i = 0; i < num_surface_elements; ++i) {
		top_file >> n0 >> n1 >> n2 >> n3;
		elem[i].set_nodes_linear(n0, n1, n2, n3);
	}
	
	// Read in interior elements
	getline(top_file, line);
	getline(top_file, line);
	for(i = 0; i < num_interior_elements; ++i) {
		top_file >> n0 >> n1 >> n2 >> n3;
		elem[i + num_surface_elements].set_nodes_linear(n0, n1, n2, n3);
	}
	
	// Close and return pointer
	top_file.close();
	return;
}

// Find nearest node in second structure to given node
int find_nearest_node(vector3 one_node, vector3 *from_node, int num_nodes_from) {
	
	int i, nearest_node;
	double shortest_distance = INFINITY;
	vector3 distance;
	for(i = 0; i < num_nodes_from; ++i) {
		distance = from_node[i] - one_node;
	//	cout << "Distance " << distance.r << endl;
	//	cout << "Shortest distance " << shortest_distance << endl;
		if(distance.r < shortest_distance) {
			shortest_distance = distance.r;
			nearest_node = i;
			cout << "Nearest node " << nearest_node << endl;
		}
	}
	return nearest_node;
}

// Find only containing element
tet_element * find_containing_element(vector3 contained_node, int nearest_node_index, tet_element *from_top, vector3 *from_node, int num_elements_from) {
	
	int i, j, success = 0, count = 0;
	double angle;
	vector3 normal[4], face_midpoint[4], diff[4];
	
	for(i = 0; i < num_elements_from; ++i) {
		
		for(j = 0; j < 4; ++j) {
			
			// Get all normal vectors and midpoints, calc_diff
			normal[j] = from_top[i].get_face_normal(j, from_node);
			face_midpoint[j] = from_top[i].get_face_midpoint(j, from_node);
			diff[j] = face_midpoint[j] - contained_node;
			diff[j].normalise();
			angle = acos(normal[j].dot(diff[j]));
			if(angle > M_PI / 2.0 || angle < - M_PI/2.0) {
				success = 0;
				break;
			} else {
				success++;
			}
		}
		
		if(success == 4) {
			count++;
			//return &from_top[i];
		}

	}
	cout << count << endl;
	return NULL;
}

int main(int argc, char **argv) {

	// Check for the right input
	if(argc != 6) {
		cout << "Usage ./" << argv[0] << " [INPUT from .node] [INPUT from .top] [INPUT to .node] [INPUT to .top] [OUTPUT .map fname]" << endl;
		exit(0);
	}

	// Extract input
	string from_node_fname, to_node_fname, from_top_fname, to_top_fname, output_map_fname;
	from_node_fname = argv[1];
	to_node_fname = argv[3];
	from_top_fname = argv[2];
	to_top_fname = argv[4];
	output_map_fname = argv[5];

	// Create node lists and topologies
	int num_elements_from, num_elements_to, num_nodes_from, num_nodes_to;
	vector3 *from_node, *to_node;
	tet_element *from_top, *to_top;
	
	num_nodes_from = get_num_nodes_from_file(from_node_fname);
	num_nodes_to = get_num_nodes_from_file(to_node_fname);
	num_elements_from = get_num_elements_from_file(from_top_fname);
	num_elements_to = get_num_elements_from_file(to_top_fname);
	
	// Extra work node in from_node
	from_node = new vector3[num_nodes_from + 1];
	to_node = new vector3[num_nodes_to];
	from_top = new tet_element[num_elements_from];
	to_top = new tet_element[num_elements_to];
	
	extract_and_create_nodes(from_node_fname, from_node);
	extract_and_create_nodes(to_node_fname, to_node);
	extract_and_create_topologies(from_top_fname, from_top);
	extract_and_create_topologies(to_top_fname, to_top);
	
	// Find nearest node, and therefore nearest element / containing element
	int i, j, nearest_node, elem_nearest_node;
	tet_element *containing_element, new_containing_element;
	vector3 translation;
	for(i = 0; i < num_nodes_to; ++i) {
		nearest_node = find_nearest_node(to_node[i], from_node, num_nodes_from);
		cout << "Nearest node = " << nearest_node << endl;
		
		// For all elements containing nearest node, see if node is in it!
		containing_element = NULL;
		containing_element = find_containing_element(to_node[i], nearest_node, from_top, from_node, num_elements_from);
		//cout << containing_element->n[0] << " " << containing_element->n[1] << " " << containing_element->n[2] << " " << containing_element->n[3] << endl;
		//exit(0);
	/*	if(containing_element == NULL) {
			
			// Get any element containing nearest node and make new element containing to_node[i] using work node
			containing_element = find_element_containing(nearest_node, from_top);
			new_containing_element.set_nodes(containing_element.n[0], containing_element.n[1], containing_element.n[2], containing_element.n[3]);
			for(j = 0; j < 4; ++j) {
				if(new_containing_element.n[j] == nearest_node) {
					translation.set_pos(to_node[i].x - from_node[nearest_node].x, to_node[i].y - from_node[nearest_node].y, to_node[i].z - from_node[nearest_node].z);7
					translation *= 2.0;
					from_node[num_nodes_from].set_pos(from_node[nearest_node].x + translation.x, from_node[nearest_node].y + translation.y, from_node[nearest_node].z + translation.z);
					new_containing_element.n[j] = num_nodes_from;
					elem_nearest_node = j;
				}
			}
			
			// Get position of node using nodes in this element
			
		} else {
			
			// Get position of node using nodes in this element
			 
		}*/
	}
	return 0;
}
