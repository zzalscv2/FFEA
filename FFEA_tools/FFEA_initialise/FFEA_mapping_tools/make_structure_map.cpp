#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include "Vectors.hpp"
#include "Elements.hpp"
#include "Matrices.hpp"

using namespace std;

enum runType {fromPDB, fromNode};

string getFileExt(const string& s) {

   size_t i = s.rfind('.', s.length());
   if (i != string::npos) {
      return(s.substr(i+1, s.length() - i));
   }

   return("");
}

int extract_nodes_from_node(string fname, vector3 *&node) {

	// Will extract nodes and move to centroid whilst it does it!
	int i, num_nodes = 0, num_interior_nodes = 0, num_surface_nodes = 0;
	double x, y, z;
	vector3 centroid;
	string line;
	
	// Open file
	ifstream nfile;
	nfile.open(fname.c_str());

	// Check correct file type
	getline(nfile, line);
	if (line.compare("ffea node file") != 0) {
		cout << "Error. Expected 'ffea node file', got '" << line << "'. Possibly not an ffea node file." << endl;
		return -1;
	}
	
	// Get num_nodes and build node_list
	nfile >> line >> num_nodes;
	cout << "\tFound " << num_nodes << " nodes." << endl;
	nfile >> line >> num_surface_nodes;
	nfile >> line >> num_interior_nodes;

	// Get some memory!
	node = new vector3[num_nodes];
	
	// Read in surface nodes
	getline(nfile, line);
	getline(nfile, line);
	for(i = 0; i < num_surface_nodes; ++i) {
		nfile >> x >> y >> z;
		node[i].set_pos(x, y, z);
		centroid.x += x;
		centroid.y += y;
		centroid.z += z;
	}
	
	// Read in interior nodes
	getline(nfile, line);
	getline(nfile, line);
	for(i = 0; i < num_interior_nodes; ++i) {
		nfile >> x >> y >> z;
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
	nfile.close();
	return num_nodes;
}

int extract_nodes_and_topology_from_pdb(string fname, vector3 *&node, string *&top) {

	// Will extract nodes and move to centroid whilst it does it!
	int i, num_nodes = 0;
	double x, y, z;
	vector3 centroid;
	char buf[100];
	FILE *fin;
	
	// Get number of nodes first so memory can be block allocated
	
	// Open file
	fin = fopen(fname.c_str(), "r");
	while(true) {
		if(feof(fin)) {
			break;	
		} else {
			fgets(buf,5,fin);
			if(strcmp(buf, "ATOM") == 0) {
				num_nodes++;
			}
			fgets(buf,96,fin);
		}
	}

	// Allocate some memory
	node = new vector3[num_nodes];
	top = new string[num_nodes];
	// Go to start of file
	fseek(fin, 0, SEEK_SET);

	i = -1;
	while(true) {
		if(feof(fin)) {
			break;	
		} else {
			fgets(buf,5,fin);
			if(strcmp(buf, "END") == 0) {
				break;
			} else if(strcmp(buf, "ATOM") == 0) {
				i++;
				
				// First, empty space
				fgets(buf, 10, fin);
				
				// Now, atom type
				fgets(buf, 4, fin);
				
				top[i] = string(buf);
				
				// Empty space
				fgets(buf,15,fin);
				
				// Atom pos
				fgets(buf,9,fin);
				x = atof(buf);

				centroid.x += x;
				fgets(buf,9,fin);
				y = atof(buf);
				centroid.y += y;
				fgets(buf,9,fin);
				z = atof(buf);
				centroid.z += z;

				fgets(buf,100,fin);
				node[i].set_pos(x, y, z);
			} else {
				fgets(buf,100,fin);
			}
		}
	}
	
	cout << "\tFound " << num_nodes << " atoms." << endl;
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
	fclose(fin);
	return num_nodes;
}


int extract_topology_from_top(string fname, vector3 *node, tetra_element *&elem) {

	int i, num_elements, num_surface_elements, num_interior_elements;
	int n[10];
	
	string line;
	
	// Open file
	ifstream tfile;
	tfile.open(fname.c_str());
	
	// Check correct file type
	getline(tfile, line);
	if (line.compare("ffea topology file") != 0) {
		cout << "Error. Expected 'ffea topology file', got '" << line << "'. Possibly not an ffea topology file." << endl;
		exit(0);
	}
	
	// Get num_elements and build node_list
	tfile >> line >> num_elements;
	cout << "\tFound " << num_elements << " elements." << endl;
	tfile >> line >> num_surface_elements;
	tfile >> line >> num_interior_elements;
	
	// Allocate some memory
	elem = new tetra_element[num_elements];
	
	// Read in surface elements
	getline(tfile, line);
	getline(tfile, line);
	for(i = 0; i < num_surface_elements; ++i) {
		tfile >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5] >> n[6] >> n[7] >> n[8] >> n[9];
		elem[i].set_structure(n[0], n[1], n[2], n[3], node);
	}

	// Read in interior elements
	getline(tfile, line);
	getline(tfile, line);
	for(i = 0; i < num_interior_elements; ++i) {
		tfile >> n[0] >> n[1] >> n[2] >> n[3] >> n[4] >> n[5] >> n[6] >> n[7] >> n[8] >> n[9];
		elem[i + num_surface_elements].set_structure(n[0], n[1], n[2], n[3], node);
	}

	// Close and return pointer
	tfile.close();
	return num_elements;
}

tetra_element * get_containing_element(vector3 tnode, tetra_element * belem, int num_elements) {
	
	// Check all elements until one is found which contains the node
	bool success;
	int i, j;
	vector3 nf, fc, fn, c;
	for(i = 0; i < num_elements; ++i) {
	
		// To know if inside, we must check all face centroids positively project onto the node-face vector
		// There may be a simpler way, but I'm a little lazy. Here, have a hug <======(^-^)======>
		
		// Calculate face centroids, face normals, and node-face vectors. 
		// Dot together. If negative, break. If all positive, return the element
		success = true;
		for(j = 0; j < 4; ++j) {
			fn = belem[i].get_face_normal(j, false);
			fc = belem[i].get_face_centroid(j);
			nf = fc - tnode;
			if(fn.dot(nf) < 0) {
				success = false;
				break;
			}
		}
		
		if(success) {
			return &belem[i];
		}
	}
	return NULL;
}

tetra_element * get_nearest_element(vector3 tnode, tetra_element * belem, int num_elements) {

	// Ok, so we're just going for the closest element now. Use element centroid - tnode magnitude
	int i;
	double distance = INFINITY;
	vector3 test;
	tetra_element *celem;
	for(i = 0; i < num_elements; ++i) {
		test = belem[i].centroid - tnode;
		if(test.r < distance) {
			distance = test.r;
			celem = &belem[i];
		}
	}
	return celem;
	 
}

void map_node_using_element(vector3 node, tetra_element *elem, double *map) {
	
	// We need to solve some linear equations
	int i;
	vector3 local_coord, local_coeff;
	matrix33 edgemat;
	
	local_coord = node - *(elem->node[0]);
	for(i = 1; i < 4; ++i) {
		edgemat.val[0][i - 1] = elem->node[i]->x - elem->node[0]->x;
		edgemat.val[1][i - 1] = elem->node[i]->y - elem->node[0]->y;
		edgemat.val[2][i - 1] = elem->node[i]->z - elem->node[0]->z;
	}
	
	local_coeff = edgemat.get_inverse().apply(local_coord);
	map[0] = 1;
	map[1] = local_coeff.x;
	map[2] = local_coeff.y;
	map[3] = local_coeff.z;
	for(i = 0; i < 3; ++i) {
		map[0] -= map[i + 1];
	}
}

void map_node_using_closest_nodes(vector3 node, vector<int> node_list, vector3 *bnode, double *map) {
	
	// Map should be the same size as node_list (done outside this function)
	
	// Again, need to solve some linear equations
	int i, j, factor;
	vector3 centroid, local_coord, local_coeff;
	matrix33 edgemat;
	
	// We need to know how many nodes are in each set
	factor = (node_list.size() - 1) / 3;			//elem->node[i] -> bnode[node_list.at(i)]
	local_coord = node - bnode[node_list.at(0)];

	// Build matrix using method of choice
	/*for(i = 0; i < 3; ++i) {
		for(j = 0; j < factor; ++j) {
			edgemat.val[0][i] += bnode[node_list.at((factor * i) + j + 1)].x - bnode[node_list.at(0)].x;
			edgemat.val[1][i] += bnode[node_list.at((factor * i) + j + 1)].y - bnode[node_list.at(0)].y;
			edgemat.val[2][i] += bnode[node_list.at((factor * i) + j + 1)].z - bnode[node_list.at(0)].z;
		}
	}*/
	
	for(i = 0; i < 3; ++i) {
		for(j = 0; j < factor; ++j) {
			edgemat.val[0][i]= 0.0;
			edgemat.val[1][i]= 0.0;
			edgemat.val[2][i]= 0.0;
		}
		
		for(j = 0; j < factor; ++j) {
			edgemat.val[0][i] += bnode[node_list.at(factor * i + j + 1)].x - bnode[node_list.at(0)].x;
			edgemat.val[1][i] += bnode[node_list.at(factor * i + j + 1)].y - bnode[node_list.at(0)].y;
			edgemat.val[2][i] += bnode[node_list.at(factor * i + j + 1)].z - bnode[node_list.at(0)].z;
		}
	}
	
	// Get coefficients
	local_coeff = edgemat.get_inverse().apply(local_coord);
	if(fabs(local_coeff.x) > 5 || fabs(local_coeff.y) > 5 || fabs(local_coeff.z) > 5) {
		cout << local_coeff.x << " " << local_coeff.y << " " << local_coeff.z << "        " << local_coord.get_mag() << "      " << node_list.size() << endl;
	}
	
	// Create local map
	/*map[0] = 1;
	for(i = 0; i < 3; ++i) {
		for(j = 0; j < factor; ++j) {
			if(i == 0) {
				map[factor * i + j + 1] = local_coeff.x;
				map[0] -= local_coeff.x;
			} else if (i == 1) {
				map[factor * i + j + 1] = local_coeff.y;
				map[0] -= local_coeff.y;
			} else {
				map[factor * i + j + 1] = local_coeff.z;
				map[0] -= local_coeff.z;
			}
		}
	}*/
	for(i = 0; i < node_list.size(); ++i) {
		map[i] = 1.0 / node_list.size();
	}
}

void add_local_map_to_global(double *node_map, double **global_map, int target_index, vector<int> node_list) {
	
	int i;
	for(i = 0; i < node_list.size(); ++i) {
		global_map[target_index][node_list.at(i)] = node_map[i];
	}
}

vector<int> get_nodes_within_range(vector3 tnode, vector3 *bnode, string *btop, int num_bnodes, double range) {
	
	int i, j, itemp;
	double dtemp;
	vector3 separation;
	vector<int> within_range;
	vector<double> distance;
	for(i = 0; i < num_bnodes; ++i) {
		
		// First of all, if not CA, continue
		if(btop[i] != "CA ") {
			continue;
		}
		
		separation = tnode - bnode[i];

		if(separation.r < range) {
			within_range.push_back(i);
			distance.push_back(separation.r);
			
			// Sort as we go
			for(j = distance.size() - 1; j > 0; --j) {
				if(distance.at(j) < distance.at(j - 1)) {
					dtemp = distance.at(j);
					distance.at(j) = distance.at(j - 1);
					distance.at(j - 1) = dtemp;
					
					itemp = within_range.at(j);
					within_range.at(j) = within_range.at(j - 1);
					within_range.at(j - 1) = itemp;
				}
			}
		}
	}
	//for(i = 0; i < within_range.size(); ++i) {
		//cout << within_range.at(0) << " " << distance.at(0) << endl;
	//}
	//exit(0);
	return within_range;
}

void write_map_to_file(string map_fname, double **map, int num_bnodes, int num_tnodes) {

	int i, j;

	// Open file
	ofstream map_file;
	map_file.open(map_fname.c_str());

	// Write initial stuff (always dense. external script will sparse it up)
	map_file << "FFEA Kinetic Conformation Mapping File ";
	map_file << "(Dense)" << endl;

	// Get num entries
	int num_entries = 0;
	for(i = 0; i < num_tnodes; ++i) {
		for(j = 0; j < num_bnodes; ++j) {

			// If irrelevant, don't include (maybe a better condition is how many non-zeroes on each line)
			if (fabs(map[i][j]) >= 1e-5) {
				num_entries++;
			}
		}
	}
	map_file << "num_nodes_from " << num_bnodes << endl;
	map_file << "num_nodes_to " << num_tnodes << endl;
	map_file << "num_entries " << num_entries << endl;
	map_file << "map:" << endl;
	
	for(i = 0; i < num_tnodes; ++i) {
		cout << "\r" << ((100 * i) / num_tnodes) << "\% of map written...";
		for(j = 0; j < num_bnodes; ++j) {

			// If irrelevant, don't include
			if (fabs(map[i][j]) < 1e-5) {
				map_file << 0 << " ";
			} else {
				map_file << map[i][j] << " ";
			}
		}
		map_file << endl;
	}
	cout << "\r" << 100 << "\% of map written...";
	cout << "done!" << endl;
}

int main(int argc, char **argv) {

	// Get args
	int i, j, k;
	string flags[] = {"-i", "-t", "-o", "-m", "-s", "-r"};
	string infname = "", topfname = "", targetfname = "", mapfname = "";
	double scale = 1.0, range = 10;
	double **map;
	
	// Make sure args come in pairs
	if((argc - 1) % 2 != 0 || argc == 1) {
		cout << "Usage: make_structure_map -i <Base .node/.pdb> -t <Base topology .top> -o <Target .node/.pdb> -m <Map fname (.map)> -s <scale>" << endl;
		exit(0);
	}
	
	for(i = 1; i < argc; i += 2) {
		if (flags[0].compare(argv[i]) == 0) {
			infname = argv[i + 1];
		} else if (flags[1].compare(argv[i]) == 0){
			topfname = argv[i + 1];
		} else if (flags[2].compare(argv[i]) == 0){
			targetfname = argv[i + 1];
		} else if (flags[3].compare(argv[i]) == 0){
			mapfname = argv[i + 1];
		} else if (flags[4].compare(argv[i]) == 0){
			scale = atof(argv[i + 1]);
		} else if (flags[5].compare(argv[i]) == 0){
			range = atof(argv[i + 1]);
		} else {
		
			// Ignore
			continue;
		}
	}
	//cout << infname << " " << topfname << " " << targetfname << " " << mapfname << endl;
	//exit(0);
	// Test base input files to decide how program should proceed (if .pdb / .node etc)
	string ext[] = {"node", "pdb"};
	runType prog;
	if(ext[0].compare(getFileExt(infname)) == 0) {
		prog = fromNode;
	} else if (ext[1].compare(getFileExt(infname)) == 0) {
		prog = fromPDB;
	} else {
		cout << "Error. Unrecognised extension on the input file. Cannot map from structure of type '" << getFileExt(infname) << "'" << endl;
	}
	
	// Now, test for topology if node input
	if(prog == fromNode) {
		if(topfname.compare("") == 0) {
			cout << "Error. If mapping from a FFEA node file (.node), we need a topology file (.top)." << endl;
			exit(0);
		}
	}
	
	// We're done testing the initialisation. Let's make some maps
	
	// Firstly, build target structures
	int num_tnodes = 0;
	vector3 *tnode = NULL;
	string *ttop = NULL;
	cout << "Building target structure..." << endl;
	if(ext[0].compare(getFileExt(targetfname)) == 0) {
		num_tnodes = extract_nodes_from_node(targetfname, tnode);
	} else if (ext[1].compare(getFileExt(targetfname)) == 0) {
		num_tnodes = extract_nodes_and_topology_from_pdb(targetfname, tnode, ttop);
	} else {
		cout << "Error. Unrecognised extension on the input file. Cannot map from structure of type '" << getFileExt(infname) << "'" << endl;
	}
	cout << "done!" << endl;
	
	//
	// Now, build input structures
	//
	int num_bnodes = 0, num_belements = 0;
	vector3 *bnode = NULL;
	string *btop = NULL;
	tetra_element *belem = NULL;
	
	cout << "Building base structure..." << endl;
	if(prog == fromNode) {
		num_bnodes = extract_nodes_from_node(infname, bnode);
		num_belements = extract_topology_from_top(topfname, bnode, belem);
	} else {
		num_bnodes = extract_nodes_and_topology_from_pdb(infname, bnode, btop);
		
	}
	cout << "done!" << endl;

	// Build the map object
	map = new double*[num_tnodes];
	cout << "Building map object..." << endl << flush;
	for(i = 0; i < num_tnodes; ++i) {
		cout << "\r" << ((100 * i) / num_tnodes) << "\% of map initialised...";
		map[i] = new double[num_bnodes];
		for(j = 0; j < num_bnodes; ++j) {
			map[i][j] = 0.0;
		}
	}
	cout << "\r" << 100 << "\% of map initialised...";
	cout << "done!" << endl;
	
	//
	// Mapping algorithm
	//
	
	cout << "Calculating map..." << endl;
	vector<int> node_list;
	if(prog == fromNode) {
	
		// Use existing topology to build the map
		
		// For every target node
		tetra_element *celem;
		double node_map[4]; 
		for(i = 0; i < num_tnodes; ++i) {
		
			cout << "\r" << ((100 * i) / num_tnodes) << "\% of map calculated..." << flush;
			
			// Find which element it is in, if any
			celem = get_containing_element(tnode[i], belem, num_belements);
			if(celem == NULL) {

				// Find closest element if outside of structure
				celem = get_nearest_element(tnode[i], belem, num_belements);
			}
				
			// We have an element! Work out how node is positioned wrt this element
			map_node_using_element(tnode[i], celem, node_map);
			node_list.clear();
			for(j = 0; j < 4; ++j) {
				node_list.push_back(celem->n[j]);
			}
			add_local_map_to_global(node_map, map, i, node_list);
		}

	} else {
	
		// Use closest N nodes, put them in order, create 3 sets out of them to reduce degrees of freedom, and go
		// Cutoff distance = range
		int tries;

		double *node_map;
		
		for(i = 0; i < num_tnodes; ++i) {
			
			node_list.clear();
			cout << "\r" << ((100 * i) / num_tnodes) << "\% of map calculated...";
			
			tries = 0;
			while (node_list.size() < 3) {
				node_list = get_nodes_within_range(tnode[i], bnode, btop, num_bnodes, range * (++tries));
			}
			
			// Bring down to nearest factor of 3 (using closest node as our reference point)
			while((node_list.size() - 1) % 3 != 0) {
				node_list.pop_back();
			}
			
			// Now, make a map fro these nodes
			node_map = new double[node_list.size()];
			map_node_using_closest_nodes(tnode[i], node_list, bnode, node_map);
		//	for(j = 0; j < node_list.size(); ++j) {
		//		cout << node_map[j] << endl;
		//	}
			add_local_map_to_global(node_map, map, i, node_list);
			
			delete[] node_map;
			node_map = NULL;
		}
	}
	
	// Now, write to file
	cout << "\r" << 100 << "\% of map calculated...";
	cout << "done!" << endl;

	// Print out the whole map
	cout << "Writing map to " << mapfname << "..." << endl;
	write_map_to_file(mapfname, map, num_bnodes, num_tnodes);
	cout << "done! " << endl;
	return 0;
}
