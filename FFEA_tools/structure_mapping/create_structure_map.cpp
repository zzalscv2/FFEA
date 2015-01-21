#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <list>

class vector3 {
	
	public:
		
		vector3(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
			calc_mag();
		}
		
		vector3() {
			this->x = 0.0;
			this->y = 0.0;
			this->z = 0.0;
			this->mag = 0.0;
		}
		
		~vector3() {
			this->x = 0.0;
			this->y = 0.0;
			this->z = 0.0;
			this->mag = 0.0;
		}
		
		// Overloaded operators
		vector3 operator+(vector3 b) {
			return vector3(x + b.x, y + b.y, z + b.z);
		}
		
		vector3 operator-(vector3 b) {
			return vector3(x - b.x, y - b.y, z - b.z);
		}
		
		vector3 operator*(double b) {
			return vector3(x * b, y * b, z * b);
		}
		void operator+=(vector3 b) {
			x += b.x;
			y += b.y;
			z += b.z;
		}
		
		void operator*=(double b) {
			x *= b;
			y *= b;
			z *= b;
		}
		
		void set_pos(double x, double y, double z) {
			this->x = x;
			this->y = y;
			this->z = z;
		}
		
		static vector3 cross_product(vector3 a, vector3 b) {
			vector3 c;
			c.x = a.y * b.z - a.z * b.y;
			c.y = b.x * a.z - a.x * b.z;
			c.z = a.x * b.y - a.y * b.x;
			return c;
		}
		
		static double dot_product(vector3 a, vector3 b) {
			return a.x * b.x + a.y * b.y + a.z * b.z;
		}
		
		void calc_mag() {
			mag = sqrt(x * x + y * y + z * z);
		}
		
		void normalise() {
			calc_mag();
			if(mag == 0.0) {
				return;
			}
			x /= mag;
			y /= mag;
			z /= mag;
		}
		
		// Variables
		double x, y, z, mag;
};

class matrix33{
	
	public:
		
		matrix33() {
			el = new double*[3];
			for(int i = 0; i < 3; ++i) {
				el[i] = new double[3];
			}
		}
		
		matrix33(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
			el = new double*[3];
			for(int i = 0; i < 3; ++i) {
				el[i] = new double[3];
				if(i = 0) {
					el[i][0] = xx;
					el[i][1] = xy;
					el[i][2] = xz;
				} else if(i = 1) {
					el[i][0] = yx;
					el[i][1] = yy;
					el[i][2] = yz;
				} else {
					el[i][0] = zx;
					el[i][1] = zy;
					el[i][2] = zz;
				}
			}
		}
		
		~matrix33() {
			delete[] el;
			el = NULL;
		}
		
		void apply(vector3 *a) {
			vector3 b;
			b.x = el[0][0] * a->x + el[0][1] * a->y + el[0][2] * a->z;
			b.y = el[1][0] * a->x + el[1][1] * a->y + el[1][2] * a->z;
			b.z = el[2][0] * a->x + el[2][1] * a->y + el[2][2] * a->z;
			a = &b;
		}
		
		void create_rotation_matrix(vector3 a, double angle) {
			double c = cos(angle);
			double s = sin(angle);
			
			el[0][0] = c + a.x * a.x * (1 - c);
			el[0][1] = a.x * a.y * (1 - c) - a.z * s;
			el[0][2] = a.x* a.z * (1 - c) + a.y * s;
			el[1][0] = a.y * a.x * (1 - c) + a.z * s;
			el[1][1] = c + a.y * a.y * (1 - c);
			el[1][2] = a.y * a.z * (1 - c) - a.x * s;
			el[2][0] = a.z * a.x * (1 - c) - a.y * s;
			el[2][1] = a.z * a.y * (1 - c) + a.x* s;
			el[2][2] = c + a.z * a.z * (1 - c); 	
		}
		
		void print_to_screen() {
			printf("Matrix is:\n");
			for(int i = 0; i < 3; ++i) {
				for(int j = 0; j < 3; ++j) {
					printf("%e ", el[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	private:
		
		double **el;
};
class Tetrahedron{
	
	public:
		
		Tetrahedron(int n0, int n1, int n2, int n3) {
			n_index[0] = n0;
			n_index[1] = n1;
			n_index[2] = n2;
			n_index[3] = n3;
		}
		
		// Variables
		int n_index[4];
};

class Volume {

	public:
		
		Volume() {
			num_nodes = 0;
			node = NULL;
			num_elements = 0;
			element.clear();
			centroid.set_pos(0.0, 0.0, 0.0);
		}
		
		~Volume() {
			num_nodes = 0;
			delete[] node;
			node = NULL;
			num_elements = 0;
			element.clear();
			centroid.set_pos(0.0, 0.0, 0.0);
		}
		
		int init(char *vol_fname) {
			
			int i, n0, n1, n2, n3;
			char line[255];
			FILE *vol_file;
			vol_file = fopen(vol_fname, "r");
			
			// Find element data
			for(;;) {
				fgets(line, 255, vol_file);
				if(strcmp(line, "volumeelements\n") == 0) {
					break;
				}
			}
			
			// Get num_elements
			fscanf(vol_file, "%d\n", &num_elements);

			// Get all element data
			for(i = 0; i < num_elements; ++i) {
				fscanf(vol_file, "%*d %*d %d %d %d %d\n", &n0, &n1, &n2, &n3);
				element.push_back(new Tetrahedron(n0, n1, n2, n3));
			}
			
			// Find node data
			for(;;) {
				fgets(line, 255, vol_file);
				if(strcmp(line, "points\n") == 0) {
					break;
				}
			}
			
			// Get num_nodes
			fscanf(vol_file, "%d\n", &num_nodes);
			
			// Create array of nodes
			node = new vector3[num_nodes];
			
			// Get all node data
			for(i = 0; i < num_nodes; ++i) {
				fscanf(vol_file, "%lf %lf %lf\n", &node[i].x, &node[i].y, &node[i].z);
			}
			
			fclose(vol_file);
			return 0;
			
		}
		
		vector3 get_node(int index) {
			return node[index];
		}
		
		void translate(vector3 trans_vec) {
			for(int i = 0; i < num_nodes; ++i) {
				node[i] += trans_vec;
			}
		}
		
		vector3 calc_centroid() {
			centroid.set_pos(0.0, 0.0, 0.0);
			for(int i = 0; i < num_nodes; ++i) {
				centroid += node[i];
			}
			centroid *= 1.0 / num_nodes;
			return centroid;
		}
		
		void apply_rotation(matrix33 a) {
			for(int i = 0; i < num_nodes; ++i) {
				printf("%d %lf %lf %lf\n", i, node[i].x, node[i].y, node[i].z);
				a.apply(&node[i]);
				printf("%d %lf %lf %lf\n", i, node[i].x, node[i].y, node[i].z);
				exit(0);
			}
		}
	private:
		
		std::list<Tetrahedron*> element;
		int num_elements;
		vector3 *node;
		int num_nodes;
		vector3 centroid;
		 
};

// Global function for calculating overlap of differing meshes
double calc_vol_overlap(Volume *vol1, Volume *vol2) {
	double overlap = 0.001;
	return overlap;
}

int main(int argc, char **argv) {
	
	// Import arguments
	if(argc != 5) {
		printf("Usage: %s [Stationary .vol file] [control_node_1] [.vol file to move] [move_control]\n", argv[0]);
		return -1;
	}

	// Import both sets of structures
	Volume *stat_vol = new Volume();
	stat_vol->init(argv[1]);
	int stat_control_index = atoi(argv[2]) - 1;
	Volume *vol_to_move = new Volume();
	vol_to_move->init(argv[3]);
	int move_control_index = atoi(argv[4]) - 1;

	// Move control node of two onto one
	vector3 stat_control = stat_vol->get_node(stat_control_index);
	vector3 move_control = vol_to_move->get_node(move_control_index);
	vector3 control_node_translate = stat_control - move_control;
	vol_to_move->translate(control_node_translate);
	
	// Move both to the origin!
	vector3 origin_translate = move_control * -1;
	stat_vol->translate(origin_translate);
	vol_to_move->translate(origin_translate);
	stat_control = stat_vol->get_node(stat_control_index);
	
	// Calculate centroid of both, then create vectors from control node to centroid, then calculate normal 
	vector3 stat_centroid = stat_vol->calc_centroid();
	vector3 move_centroid = vol_to_move->calc_centroid();
	vector3 stat_axial = stat_centroid - stat_control;
	vector3 move_axial = move_centroid - move_control;
	vector3 rotation_vector = vector3::cross_product(stat_centroid, move_centroid);
	rotation_vector.normalise();
	
	// Rotate about normal. Iterate to maximum volume overlap
	double new_overlap = INFINITY, old_overlap = 0.0, test_overlap, rot_angle, current_angle = 0.0, ang_lower_lim = 0.0, ang_upper_lim = 2 * M_PI, highest_overlap = 0.0;
	matrix33 rot_mat;
	int runs = -1;
	double threshold = 0.01;
	while(runs++ < 100) {
		// Check if overlap is satisfactory (on first run it won't be!)
		if((new_overlap - old_overlap) / new_overlap < threshold && runs != 0) {
			break;
		}
		
		// If not, iterate again by splitting current sector into 10 parts and checking the overlap at each one.
		// We want the highest overlap, and our new section would be between the two either side of it
		highest_overlap = 0.0;
		rot_angle = ang_upper_lim - ang_lower_lim / 10.0;
		rot_mat.create_rotation_matrix(rotation_vector, rot_angle);
		for(int i = 0; i < 10; ++i) {
			printf("Hi1");
			vol_to_move->apply_rotation(rot_mat);
			printf("Hi2");
			current_angle += rot_angle;
			test_overlap = calc_vol_overlap(stat_vol, vol_to_move);
			printf("Hi3");
			if(test_overlap > highest_overlap) {
				highest_overlap = test_overlap;
				ang_lower_lim = current_angle - rot_angle;
				ang_upper_lim = current_angle + rot_angle;
			}
			printf("Hi4");
		}
		
		// Rotate back to low_limit ovelap position
		rot_mat.create_rotation_matrix(rotation_vector, ang_lower_lim - current_angle);
		vol_to_move->apply_rotation(rot_mat);
		current_angle = ang_lower_lim;
		old_overlap = new_overlap;
		new_overlap = test_overlap;
	}
	
	// Rotate to the high_overlap position
	rot_mat.create_rotation_matrix(rotation_vector, ((ang_upper_lim + ang_lower_lim) / 2.0) - current_angle);
	vol_to_move->apply_rotation(rot_mat);
	current_angle = (ang_upper_lim + ang_lower_lim) / 2.0;
	
	// Calculate axial vector form control node to new centroid
	//stat_centroid = 
	// Rotate about axial vector. Iterate to maximum volume overlap
	
	// Create map
	
	// Print shit out
}
