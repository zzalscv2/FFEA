#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <vector>

#define ERROR   "Error reading line"

#define MAX_FNAME_LENGTH 255

#define MAP_INDEX(X, Y, Z) (((Z)*map->NY*map->NX) + ((Y)*map->NX) + ((X)))

// As defined by the spec
#define HEADER_SIZE_IN_WORDS 256
#define WORD_SIZE 4

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

int read_pdb_frame(FILE *pdb_file_stream, int num_atoms, vector3 *frame)
{
	const int max_line_length = 100;
	char line[max_line_length];
	char record_name[4];

	int i = 0;
	float x, y, z;
        for(;;) {
		if(fgets(line, max_line_length, pdb_file_stream) != NULL) {
			// get first 4 characters of line
			for(int j = 0; j < 4; j++) {
				record_name[j] = line[j];
			}
			if(strstr(record_name, "ATOM") != NULL) {
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
			} 
		} else {
			printf("End of .pdb reached\n");
			return i;
		}
	}
}

int get_num_atoms(char *pdb_fname)
{
	const int max_line_length = 100;
	char line[max_line_length];
	char record_name[4];

	FILE *pdb_file_stream = NULL;
	if((pdb_file_stream = fopen(pdb_fname, "r")) == NULL) {
		printf("Could not open pdb file '%s'\n", pdb_fname);
		return -1;
	}

	int i = 0;
        for(;;) {
		if(fgets(line, max_line_length, pdb_file_stream) != NULL) {
			// get first 4 characters of line
			for(int j = 0; j < 4; j++) {
				record_name[j] = line[j];
			}
			if(strstr(record_name, "ATOM") != NULL) {
				i++;
				for(;;) {
					if(fgets(line, max_line_length, pdb_file_stream) == NULL) {
						return i;
					}
					for(int j = 0; j < 4; j++) {
						record_name[j] = line[j];
					}
					if(strstr(record_name, "ATOM") != NULL) {
						i++;
					}
				}
			}
		} else {
			return 0;
		}
	}
}

float sphere_density_func(float dist, float radius)
{
	if(dist > radius) {
		return 0;
	} else {
		return (1 - dist/radius);
	}
}

void add_density_sphere_to_map(float radius, float x, float y, float z, map_data *map)
{
	double temp = 0.0;

	// find the size of the box (in number of voxels) that contains the sphere
	int rad_vox_x_plus = ceil(radius/map->voxel_sizeX), rad_vox_y_plus = ceil(radius/map->voxel_sizeY), rad_vox_z_plus = ceil(radius/map->voxel_sizeZ);
	int rad_vox_x_minus = rad_vox_x_plus, rad_vox_y_minus = rad_vox_y_plus, rad_vox_z_minus = rad_vox_z_plus;

	// find the cell this sphere is centred on
	int vox_x = floor(x/map->voxel_sizeX), vox_y = floor(y/map->voxel_sizeY), vox_z = floor(z/map->voxel_sizeZ);

	// Define some local co-ordinates
	double local_x = x - vox_x * map->voxel_sizeX, local_y = y - vox_y * map->voxel_sizeY, local_z = z - vox_z * map->voxel_sizeZ;

	// recalculate rad_vox_* due to position within voxel
	temp = map->voxel_sizeX * (1 - rad_vox_x_minus) + radius;
	if(local_x > temp) {
		rad_vox_x_minus--;
	} else if(local_x > map->voxel_sizeX - temp) {
		rad_vox_x_plus++;
	}

	temp = map->voxel_sizeY * (1 - rad_vox_y_minus) + radius;
	if(local_y > temp) {
		rad_vox_y_minus--;
	} else if(local_y > map->voxel_sizeY - temp) {
		rad_vox_y_plus++;
	}

	temp = map->voxel_sizeZ * (1 - rad_vox_z_minus) + radius;
	if(local_z > temp) {
		rad_vox_z_minus--;
	} else if(local_z > map->voxel_sizeZ - temp) {
		rad_vox_z_plus++;
	}

	// iterate through each voxel affected by the sphere and calculate its density contribution to that voxel
	for(int i = vox_x - rad_vox_x_minus; i < vox_x + rad_vox_x_plus; i++) {
		for(int j = vox_y - rad_vox_y_minus; j < vox_y + rad_vox_y_plus; j++) {
			for(int k = vox_z - rad_vox_z_minus; k < vox_z + rad_vox_z_plus; k++) {
				if(i == -1) {
					i++;
				}
				if(j == -1) {
					j++;
				}
				if(k == -1) {
					k++;
				}

				float dx = x - i * map->voxel_sizeX;
				float dy = y - j * map->voxel_sizeY;
				float dz = z - k * map->voxel_sizeZ;
				float r = dx * dx + dy * dy + dz * dz;
				map->data[MAP_INDEX(i, j, k)] += sphere_density_func(r, radius);
				if (i < 0 || j < 0 || k < 0) {
					printf("FUCK %d %d %d\n", i, j, k);
				}
			}
		}
	}
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

	for(int i = 0; i < map->num_voxels; i++) {
		fwrite(&(map->data[i]), sizeof(float), 1, out);
	}

	return 0;
}

void get_extent(vector3 *frame, int num_atoms, float &min_x, float &min_y, float &min_z, float &max_x, float &max_y, float &max_z)
{
	min_x = INFINITY;
	min_y = INFINITY;
	min_z = INFINITY;
	max_x = -INFINITY;
	max_y = -INFINITY;
	max_z = -INFINITY;
	for(int i = 0; i < num_atoms; i++) {
		float x = frame[i].x, y = frame[i].y, z = frame[i].z;
		if(x < min_x) {
			min_x = x;
		}
		if(x > max_x) {
			max_x = x;
		}
		if(y < min_y) {
			min_y = y;
		}
		if(y > max_y) {
			max_y = y;
		}
		if(z < min_z) {
			min_z = z;
		}
		if(z > max_z) {
			max_z = z;
		}
	}
}

void translate_frame(vector3 *frame, int num_atoms, float t_x, float t_y, float t_z)
{
	for(int i = 0; i < num_atoms; i++) {
		frame[i].x += t_x;
		frame[i].y += t_y;
		frame[i].z += t_z;
	}
}

void translate_extent(float &min_x, float &min_y, float &min_z, float &max_x, float &max_y, float &max_z)
{
	max_x -= min_x;
	max_y -= min_y;
	max_z -= min_z;
	min_x = 0.0;
	min_y = 0.0;
	min_z = 0.0;
	printf("New extent: \tmin_x = %f  min_y = %f  min_z = %f\n\t\tmax_x = %f  max_y = %f max_z = %f\n", min_x, min_y, min_z, max_x, max_y, max_z);
}

int main(int argc, char **argv)
{
	char pdb_fname[MAX_FNAME_LENGTH], out_fname[MAX_FNAME_LENGTH];

	if(argc != 7) {
		printf("Usage: ./convert_pdb_to_map INPUT_PDB_FNAME OUTPUT_MAP_FNAME NX NY NZ ATOMIC_RADIUS(angstroms)\n");
		return 1;
	}

	sprintf(pdb_fname, "%s", argv[1]);
	sprintf(out_fname, "%s", argv[2]);

	int nx = atoi(argv[3]);
	int ny = atoi(argv[4]);
	int nz = atoi(argv[5]);
	float radius = atoi(argv[6]);

	// Read in first frame from pdb file
	int num_atoms = get_num_atoms(pdb_fname);

	if (num_atoms == -1) {
		printf("File not found. Please provide an appropriate pdb file.");
		return 1;

	} else if (num_atoms == 0) {
		printf("Error. Found 0 atoms in '%s'. Possibly incorrectly formatted for this tool. Please contact developers.\n", pdb_fname);
		return 1;
	}

	printf("Found %d atoms in pdb file.\n", num_atoms);

	printf("Reading first frame of pdb file...\n");
	FILE *pdb_in, *map_out;
	if((pdb_in = fopen(pdb_fname, "r")) == NULL) {
		printf("Could not open pdb file '%s'\n", pdb_fname);
		return 1;
	}

	printf("num_atoms = %d\n", num_atoms);
	vector3 *frame = new vector3[num_atoms];
	int num_atoms_read = read_pdb_frame(pdb_in, num_atoms, frame);
	printf("Read %d atoms from pdb.\n", num_atoms_read);
	fclose(pdb_in);
	printf("Done.\n");

	if((map_out = fopen(out_fname, "w")) == NULL) {
		printf("Could not open file '%s' for writing\n", out_fname);
		return 1;
	}

	// Get size of box necessary to contain protein
	float min_x, min_y, min_z, max_x, max_y, max_z;
	get_extent(frame, num_atoms, min_x, min_y, min_z, max_x, max_y, max_z);

	// Add margin of size "radius" so that edge atoms can fit in box
	min_x -= radius;
	min_y -= radius;
	min_z -= radius;
	max_x += radius;
	max_y += radius;
	max_z += radius;
	printf("Protein extent: min_x = %f min_y = %f min_z = %f\n\t\tmax_x = %f max_y = %f max_z = %f\n", min_x, min_y, min_z, max_x, max_y, max_z);

	// Box width, height, depth
	float xlength = max_x - min_x, ylength = max_y - min_y, zlength = max_z - min_z;

	// Translate atoms so that min side of box starts at (0, 0, 0)
	translate_frame(frame, num_atoms, -min_x, -min_y, -min_z);
	translate_extent(min_x, min_y, min_z, max_x, max_y, max_z);

	// Make a new map and header
	map_data map;
	map.num_voxels = nx * ny * nz;
	map.NX = nx;
	map.NY = ny;
	map.NZ = nz;
	map.voxel_sizeX = xlength/nx;
	map.voxel_sizeY = ylength/ny;
	map.voxel_sizeZ = zlength/nz;
	if(map.voxel_sizeX > 2 * radius || map.voxel_sizeY > 2 * radius || map.voxel_sizeZ > 2 * radius) {
		printf("Voxel is larger than atom. This may not end well. Change atomic radius or number of voxels maybe?\n");	
	}
	map.data = new float[map.num_voxels];

	header_info header;
	header.NC = nx;
	header.NR = ny;
	header.NS = nz;
	header.MODE = 2;
	header.NCSTART = 0;
	header.NRSTART = 0;
	header.NSSTART = 0;
	header.NX = nx;
	header.NY = ny;
	header.NZ = nz;
	header.X_LENGTH = xlength;
	header.Y_LENGTH = ylength;
	header.Z_LENGTH = zlength;
	header.ALPHA = 90.0;
	header.BETA = 90.0;
	header.GAMMA = 90.0;
	header.MAPC = 1;
	header.MAPR = 2;
	header.MAPS = 3;
	header.AMIN = 0;
	header.AMAX = 0;
	header.AMEAN = 0;
	header.ISPG = 0;
	header.NSYMBT = 0;
	header.LSKFLG = 0;
	for(int i = 0; i < 9; i++) {
		header.SKWMAT[i] = 0.0;
	}
	for(int i = 0; i < 3; i++) {
		header.SKWTRN[i] = 0.0;
	}
	for(int i = 0; i < 15; i++) {
		header.EXTRA[i] = 0.0;
	}
	header.MAP[0] = 'M';
	header.MAP[1] = 'A';
	header.MAP[2] = 'P';
	header.MAP[3] = '\0';
	header.MACHST = 0;
	header.RMS = 0;
	header.NLABL = 1;
	sprintf(header.LABEL_N, "Created by FFEA_tools: convert_pdb_to_emdb_map %s %s %d %d %d\n", pdb_fname, out_fname, nx, ny, nz);

	// Convert protein to density map
	for(int i = 0; i < num_atoms; i++) {
		add_density_sphere_to_map(radius, frame[i].x, frame[i].y, frame[i].z, &map);
	}

	// write new map to file
	write_map_to_file(&map, &header, out_fname);

	return 0;
}
