#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MersenneTwister.h"

#define ERROR	"Error reading line"

class vector3
{
	public:
		vector3()
		{
			x = 0; y = 0; z = 0;
			name = '\0';
		}
		vector3(float x, float y, float z)
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		void set(float x, float y, float z, char name)
		{
			this->x = x;
			this->y = y;
			this->z = z;
			this->name = name;
		}

		void set(vector3 *v)
		{
			this->x = v->x;
			this->y = v->y;
			this->z = v->z;
			this->name = v->name;
		}

		float x, y, z;
		char name;
};

void read_pdb_frame(FILE *pdb_file_stream, int num_atoms, vector3 *frame)
{
	const int max_line_length = 100;
	char line[max_line_length];

	int i = 0;
	float x, y, z;
        for(;;) {
		if(fgets(line, max_line_length, pdb_file_stream) != NULL) {
			if(strstr(line, "REMARK") != NULL) {
				printf("Skipping 'REMARK'\n");
				continue;
			}
			if(strstr(line, "MODEL") != NULL) {
				printf("Skipping 'MODEL'\n");
			continue;
			}
			if(strstr(line, "ATOM") != NULL) {
				printf("Reading ATOM block\n");
				for(int n = 0; n < num_atoms; n++) {
//					int check = sscanf(line, "%*s %*d %*s %*s %*d %f%f%f", &x, &y, &z);
//					int check = sscanf(line, "%*30c%8f%8f%8f", &x, &y, &z);

					char sx[9], sy[9], sz[9];
					int check = sscanf(line, "%*30c%8c%8c%8c", sx, sy, sz);
					if(check != 3) {
						printf("Error reading line in frame at atom %d\n", i);
						throw ERROR;
					}

					// work around for faulty 8 width field handling
					sx[8] = '\0';
					sy[8] = '\0';
					sz[8] = '\0';
					x = atof(sx);
					y = atof(sy);
					z = atof(sz);

					frame[i].set(x, y, z, line[13]);
					i++;
					if(fgets(line, max_line_length, pdb_file_stream) == NULL) {
						return;
					}
					if(strstr(line, "ATOM") == NULL) {
						return;
					}
				}
				return;
			}
		}
	}
}

int get_num_atoms(char *pdb_fname)
{
	const int max_line_length = 100;
	char line[max_line_length];

	FILE *pdb_file_stream = NULL;
	if((pdb_file_stream = fopen(pdb_fname, "r")) == NULL) {
		printf("Could not open pdb file '%s'\n", pdb_fname);
		return 1;
	}

	int i = 0;
        for(;;) {
		if(fgets(line, max_line_length, pdb_file_stream) != NULL) {
			if(strstr(line, "REMARK") != NULL) {
				continue;
			}
			if(strstr(line, "MODEL") != NULL) {
				continue;
			}
			if(strstr(line, "ATOM") != NULL) {
				i++;
				for(;;) {
					if(fgets(line, max_line_length, pdb_file_stream) == NULL) {
						return i;
					}
					if(strstr(line, "ATOM") != NULL) {
						i++;
					} else {
						return i;
					}
				}
			}
		} else {
			return 0;
		}
	}
}

void get_centroid(int num_atoms, vector3 *frame, vector3 *centroid)
{
	centroid->x = 0;
	centroid->y = 0;
	centroid->z = 0;
	for(int i = 0; i < num_atoms; i++) {
		centroid->x += frame[i].x;
		centroid->y += frame[i].y;
		centroid->z += frame[i].z;
	}
	centroid->x /= num_atoms;
	centroid->y /= num_atoms;
	centroid->z /= num_atoms;
}

void construct_rotation_matrix(double rotation[][3], double alpha, double beta, double gamma)
{
	rotation[0][0] = cos(alpha) * cos(beta);
	rotation[0][1] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma);
	rotation[0][2] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma);

	rotation[1][0] = sin(alpha) * cos(beta);
	rotation[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma);
	rotation[1][2] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma);

	rotation[2][0] = -sin(beta);
	rotation[2][1] = cos(beta) * sin(gamma);
	rotation[2][2] = cos(beta) * cos(gamma);
}

void mat_mult(double M[][3], vector3 *v, vector3 *result)
{
        result->x = M[0][0]*v->x + M[0][1]*v->y + M[0][2]*v->z;
        result->y = M[1][0]*v->x + M[1][1]*v->y + M[1][2]*v->z;
        result->z = M[2][0]*v->x + M[2][1]*v->y + M[2][2]*v->z;
}

void translate(vector3 *p, vector3 *translate_by)
{
	p->x += translate_by->x;
	p->y += translate_by->y;
	p->z += translate_by->z;
}

double get_diff_penalty(vector3 *frame_large, int num_atoms_frame_large, vector3 *frame_small, int num_atoms_frame_small, double dx, double dy, double dz, vector3 *centroid_large, vector3 *centroid_small, double yaw, double pitch, double roll)
{
	double rotation[3][3];
	construct_rotation_matrix(rotation, yaw, pitch, roll);
	double E = 0;
	vector3 new_p;
	vector3 trans_vec = vector3(centroid_small->x + dx, centroid_small->y + dy, centroid_small->z + dz);
	for(int p = 0; p < num_atoms_frame_small; p++) {
		mat_mult(rotation, &frame_small[p], &new_p);
		translate(&new_p, &trans_vec);
		for(int q = 0; q < num_atoms_frame_large; q++) {
			double diff_x = new_p.x - (frame_large[q].x + centroid_large->x);
			double diff_y = new_p.y - (frame_large[q].y + centroid_large->y);
			double diff_z = new_p.z - (frame_large[q].z + centroid_large->z);
			E += sqrt(diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
//			E += diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
		}
	}
	return E;
}

int main(int argc, char **argv)
{
	if(argc != 6 && argc != 9) {
		printf("./align_point_cloud [INPUT LARGE POINT CLOUD PDB] [INPUT SMALL POINT CLOUD PDB] [OUTPUT SMALL POINT CLOUD PDB] [SOLUTION FNAME] [NUM ITERATIONS] OPTIONAL{[INITIAL YAW] [INITIAL PITCH] [INITIAL ROLL]}\n");
		return 1;
	}

	int num_iter = atoi(argv[5]);
	double MAX_TRANSLATION_STEP = 5.0;
	double MAX_ROTATION_STEP = M_PI/90.0;

	MTRand rng = MTRand();

	int num_atoms_large = get_num_atoms(argv[1]);
	int num_atoms_small = get_num_atoms(argv[2]);

	printf("Num atoms in large:%d\n", num_atoms_large);
	printf("Num atoms in small:%d\n", num_atoms_small);

	vector3 frame_large[num_atoms_large];
	vector3 frame_small[num_atoms_small];

	FILE *pdb_large_file, *pdb_small_file, *pdb_out;
	if((pdb_large_file = fopen(argv[1], "r")) == NULL) {
		printf("Could not open pdb file '%s'\n", argv[1]);
		return 1;
	}

	if((pdb_small_file = fopen(argv[2], "r")) == NULL) {
		printf("Could not open pdb file '%s'\n", argv[2]);
		return 1;
	}


	printf("Reading frame from large\n");
	read_pdb_frame(pdb_large_file, num_atoms_large, frame_large);
	printf("Reading frame from small\n");
	read_pdb_frame(pdb_small_file, num_atoms_small, frame_small);

	fclose(pdb_large_file);
	fclose(pdb_small_file);

	printf("Ignoring atoms named 'N', performing fit with atoms named 'Y'...\n");
	vector3 frame_small_active[num_atoms_small];
	int num_Y_atoms = 0;
	for(int i = 0; i < num_atoms_small; i++) {
		if(frame_small[i].name == 'Y') {
			frame_small_active[num_Y_atoms].set(&frame_small[i]);
			num_Y_atoms++;
		}
	}
	printf("%d Y atoms found in small point cloud pdb.\n", num_Y_atoms);

	// use the centroid of the full cloud (since we want the rotation for the whole structure, not just the selection)
	// but the energy penalty (calculated later) is always done only with the 'active' atoms (named Y in the pdb)
	vector3 centroid_large, centroid_small;
	get_centroid(num_atoms_large, frame_large, &centroid_large);
//	get_centroid(num_atoms_small, frame_small, &centroid_small);
	get_centroid(num_Y_atoms, frame_small_active, &centroid_small);



	double dx = centroid_large.x - centroid_small.x;
	double dy = centroid_large.y - centroid_small.y;
	double dz = centroid_large.z - centroid_small.z;

	double yaw = 0;
	double pitch = 0;
	double roll = 0;

	if(argc == 9) {
		yaw = atof(argv[6]);
		pitch = atof(argv[7]);
		roll = atof(argv[8]);
	}

	vector3 trans_l = vector3(-centroid_large.x, -centroid_large.y, -centroid_large.z);
	vector3 trans_s = vector3(-centroid_small.x, -centroid_small.y, -centroid_small.z);

	// centre the clouds on the origin
	for(int i = 0; i < num_atoms_large; i++) {
		translate(&frame_large[i], &trans_l);
	}
	for(int i = 0; i < num_atoms_small; i++) {
		translate(&frame_small[i], &trans_s);
	}
	for(int i = 0; i < num_Y_atoms; i++) {
		translate(&frame_small_active[i], &trans_s);
	}

	// get starting "energy" of the configuration (using only 'active' atoms)
	double E = get_diff_penalty(frame_large, num_atoms_large, frame_small_active, num_Y_atoms, dx, dy, dz, &centroid_large, &centroid_small, yaw, pitch, roll);
	printf("Difference penalty E = %f\n", E);

	int num_improvements = 0;
	double ddx, ddy, ddz, dyaw, dpitch, droll, new_E, dE;
	for(int i = 0; i < num_iter; i++) {
		if(i%100 == 0) {
			printf("%d / %d (%d)\n", i, num_iter, num_improvements);
			printf("\033[F\033[J"); // Replace previous "progress" line

			if(num_improvements == 0) {
				MAX_TRANSLATION_STEP *= .9;
				MAX_ROTATION_STEP *= .9;
			}
			num_improvements = 0;
		}

		// Make a random translation and a random rotation
		ddx = (rng.rand() - .5) * MAX_TRANSLATION_STEP;
		ddy = (rng.rand() - .5) * MAX_TRANSLATION_STEP;
		ddz = (rng.rand() - .5) * MAX_TRANSLATION_STEP;
		dyaw = (rng.rand() - .5) * MAX_ROTATION_STEP;
		dpitch = (rng.rand() - .5) * MAX_ROTATION_STEP;
		droll = (rng.rand() - .5) * MAX_ROTATION_STEP;

		new_E = get_diff_penalty(frame_large, num_atoms_large, frame_small_active, num_Y_atoms, dx + ddx, dy + ddy, dz + ddz, &centroid_large, &centroid_small, yaw + dyaw, pitch + dpitch, roll + droll);
		dE = new_E - E;
//		p = exp(-dE/T)
//		if random.rand() < p:
		if(dE < 0) {
			num_improvements++;

			E = new_E;
			dx += ddx;
			dy += ddy;
			dz += ddz;
			yaw += dyaw;
			pitch += dpitch;
			roll += droll;
//			T -= max_T/20000;
		}
	}

	FILE *soln_out = NULL;
	if((soln_out = fopen(argv[4], "w")) == NULL) {
		printf("Could not open pdb file '%s'\n", argv[4]);
		return 1;
	}
	fprintf(soln_out, "%f %f %f %f %f %f\n", dx, dy, dz, yaw, pitch, roll);
	fclose(soln_out);

	if((pdb_out = fopen(argv[3], "w")) == NULL) {
		printf("Could not open pdb file '%s'\n", argv[3]);
		return 1;
	}

	double rotation[3][3];
	construct_rotation_matrix(rotation, yaw, pitch, roll);
	vector3 new_p;
	vector3 trans_vec = vector3(centroid_small.x + dx, centroid_small.y + dy, centroid_small.z + dz);
	for(int p = 0; p < num_atoms_small; p++) {
		mat_mult(rotation, &frame_small[p], &new_p);
		translate(&new_p, &trans_vec);
		fprintf(pdb_out, "ATOM%7d  N   GLY     1    %8.3f%8.3f%8.3f\n", p, new_p.x, new_p.y, new_p.z);
	}
	fprintf(pdb_out, "TER\nENDMDL\n");
	fclose(pdb_out);

	printf("Translation: (%f %f %f)\n", dx, dy, dz);
	printf("Rotation: (%f %f %f)\n", yaw, pitch, roll);
	printf("Rotation matrix:\n");
	for(int i = 0; i < 3; i++) {
		printf("(%f %f %f)\n", rotation[i][0], rotation[i][1], rotation[i][2]);
	}
}
