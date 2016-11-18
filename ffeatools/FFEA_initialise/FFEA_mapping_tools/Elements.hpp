#ifndef ELEMENTS_HPP
#define ELEMENTS_HPP

class tetra_element
{
	public:
		tetra_element() {
			for(int i = 0; i < 4; ++i) {
				n[i] = -1;
				node[i] = NULL;
				centroid.set_pos(0.0,0.0,0.0);
			}
		}
		
		~tetra_element() {
			for(int i = 0; i < 4; ++i) {
				n[i] = 0;
				node[i] = NULL;
				centroid.set_pos(0.0,0.0,0.0);
			}
		}
		
		void set_structure(int n0, int n1, int n2, int n3, vector3 *node_list, int ind) {
		
			// Store indices and also pointers to actual node vectors too!
			n[0] = n0;
			n[1] = n1;
			n[2] = n2;
			n[3] = n3;
			
			for(int i = 0; i < 4; ++i) {
				node[i] = &node_list[n[i]];
			}
			calc_centroid();
			index = ind;
		}
		
		void calc_centroid() {
			for(int i = 0; i < 4; ++i) {
				centroid += *node[i];
			}
			centroid *= 0.25;
		}
		
		vector3 get_face_normal(int index, bool normed) {
			vector3 normal, edge[2];
			
			// Indices: 0-012, 1-132, 2-230, 3-310
			if(index == 0) {
				edge[0] = *node[1] - *node[0];
				edge[1] = *node[2] - *node[0];
			} else if (index == 1) {
				edge[0] = *node[3] - *node[1];
				edge[1] = *node[2] - *node[1];
			} else if (index == 2) {
				edge[0] = *node[3] - *node[2];
				edge[1] = *node[0] - *node[2];
			} else if (index == 3) {
				edge[0] = *node[1] - *node[3];
				edge[1] = *node[0] - *node[3];
			}
			
			normal = edge[0].cross(edge[1]);
			//cout << normal.x << " " << normal.y << " " << normal.z << endl;
			if (normed) {
				normal *= 1.0/normal.get_mag();
			}
			
			return normal;
		}
		
		vector3 get_face_centroid(int index) {
			
			// Indices: 0-012, 1-132, 2-230, 3-310
			if(index == 0) {
				return (*node[0] + *node[1] + *node[2]) * (1.0/3.0);
			} else if (index == 1) {
				return (*node[1] + *node[3] + *node[2]) * (1.0/3.0);
			} else if (index == 2) {
				return (*node[2] + *node[3] + *node[0]) * (1.0/3.0);
			} else if (index == 3) {
				return (*node[3] + *node[1] + *node[0]) * (1.0/3.0);
			} else {
				return (*node[0]) * 0.0;
			} 
		}
		
		static bool connected(tetra_element a, tetra_element b) {
			int check = 0;
			for(int i = 0; i < 4; ++i) {
				for(int j = 0; j < 4; ++j) {
					if(a.n[i] == b.n[j]) {
						check++;
						continue;
					}
				}
			}
			if(check == 3) {
				return true;
			}
			return false;
		}

		// Data members
		int n[4];
		int index;
		vector3 *node[4];
		vector3 centroid;
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

#endif
