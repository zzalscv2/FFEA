// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <vector>

#define MAX_FNAME_LENGTH 255

#define MAP_INDEX(X, Y, Z) (((Z)*map->NY*map->NX) + ((Y)*map->NX) + ((X)))
#define MAP1_INDEX(X, Y, Z) (((Z)*map_1->NY*map_1->NX) + ((Y)*map_1->NX) + ((X)))
#define MAP2_INDEX(X, Y, Z) (((Z)*map_2->NY*map_2->NX) + ((Y)*map_2->NX) + ((X)))
#define MAP_OUT_INDEX(X, Y, Z) (((Z)*map_out->NY*map_out->NX) + ((Y)*map_out->NX) + ((X)))

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

// Extract the header and voxel data from the given emdb binary map file
int extract_data(char *map_fname, header_info *header, map_data *map)
{
	FILE *mapfile = NULL;
	int i;
	int shutup;

	if((mapfile = fopen(map_fname, "r")) == NULL) {
		printf("Could not open map file %s\n", map_fname);
		return 1;
	}

	shutup = fread(&(header)->NC, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NR, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NS, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->MODE, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NCSTART, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NRSTART, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NSSTART, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NX, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NY, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NZ, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->X_LENGTH, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->Y_LENGTH, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->Z_LENGTH, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->ALPHA,  sizeof(float), 1, mapfile);
	shutup = fread(&(header)->BETA, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->GAMMA, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->MAPC, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->MAPR, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->MAPS, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->AMIN, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->AMAX, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->AMEAN, sizeof(float), 1, mapfile);
	shutup = fread(&(header)->ISPG, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->NSYMBT, sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header)->LSKFLG, sizeof(int32_t), 1, mapfile);

	for(i = 0; i < 9; i++) {
		shutup = fread(&(header->SKWMAT[i]), sizeof(float), 1, mapfile);
	}

	for(i = 0; i < 3; i++) {
		shutup = fread(&(header->SKWTRN[i]), sizeof(float), 1, mapfile);
	}

	for(i = 0; i < 15; i++) {
		shutup = fread(&(header->EXTRA[i]), sizeof(int32_t), 1, mapfile);
	}

	header->MAP[0] = 'B'; header->MAP[1] = 'A'; header->MAP[2] = 'D'; header->MAP[3] = '\0';
	shutup = fread(header->MAP, sizeof(char), 4, mapfile);
	header->MAP[4] = '\0';

	shutup = fread(&(header->MACHST), sizeof(int32_t), 1, mapfile);
	shutup = fread(&(header->RMS), sizeof(float), 1, mapfile);
	shutup = fread(&(header->NLABL), sizeof(int32_t), 1, mapfile);

	shutup = fread(header->LABEL_N, sizeof(char), 200 * WORD_SIZE, mapfile);
	header->LABEL_N[200 * WORD_SIZE] = '\0';

//	shutup = fread(header->LABEL_N, sizeof(char), 1, mapfile);

	// Now read in the density data
/*
	map->num_voxels = header->NX * header->NY * header->NZ;
	if((header->NX != header->NY) || (header->NX != header->NZ)) {
		printf("Error: NX, NY and NZ should all be the same\n");
		return 1;
	}
	map->N = header->NX;
	map->voxel_size = header->X_LENGTH/(float)map->N;
	map->data = new float[map->num_voxels];
	shutup = fread(map->data, sizeof(float), map->num_voxels, mapfile);
*/

	if((header->NX < 1) || (header->NY < 1) || (header->NZ < 1)) {
		printf("Error: Map file has dimensions %dx%dx%d\n", header->NX, header->NY, header->NZ);
		return 1;
	}

	map->num_voxels = (header->NX + 2) * (header->NY + 2) * (header->NZ + 2);
//	if((header->NX != header->NY) || (header->NX != header->NZ)) {
//		printf("Error: NX, NY and NZ should all be the same\n");
//		return 1;
//	}
	map->NX = header->NX;
	map->NY = header->NY;
	map->NZ = header->NZ;
	map->voxel_sizeX = header->X_LENGTH/(float)(map->NX);
	map->voxel_sizeY = header->Y_LENGTH/(float)(map->NY);
	map->voxel_sizeZ = header->Z_LENGTH/(float)(map->NZ);
	map->data = new float[map->num_voxels];
	for(i = 0; i < map->num_voxels; i++) {
		map->data[i] = -INFINITY;
	}
	for(int z = 0; z < map->NZ; z++) {
		for(int y = 0; y < map->NY; y++) {
			for(i = MAP_INDEX(0,y,z); i < MAP_INDEX(header->NX,y,z); i++) {
				shutup = fread(&map->data[i], sizeof(float), 1, mapfile);
			}
		}
	}

	fclose(mapfile);

	return 0;
}

void print_header_to_file(FILE *f, header_info *header)
{
	int i;

	fprintf(f, "#\n# Header:\n");
	fprintf(f, "# NC = %d\n", header->NC);
	fprintf(f, "# NR = %d\n", header->NR);
	fprintf(f, "# NS = %d\n", header->NS);
	fprintf(f, "# MODE = %d\n", header->MODE);
	fprintf(f, "# NCSTART = %d\n", header->NCSTART);
	fprintf(f, "# NRSTART = %d\n", header->NRSTART);
	fprintf(f, "# NSSTART = %d\n", header->NSSTART);
	fprintf(f, "# NX = %d\n", header->NX);
	fprintf(f, "# NY = %d\n", header->NY);
	fprintf(f, "# NZ = %d\n", header->NZ);
	fprintf(f, "# X_LENGTH = %f\n", header->X_LENGTH);
	fprintf(f, "# Y_LENGTH = %f\n", header->Y_LENGTH);
	fprintf(f, "# Z_LENGTH = %f\n", header->Z_LENGTH);
	fprintf(f, "# ALPHA = %f\n", header->ALPHA);
	fprintf(f, "# BETA = %f\n", header->BETA);
	fprintf(f, "# GAMMA = %f\n", header->GAMMA);
	fprintf(f, "# MAPC = %d\n", header->MAPC);
	fprintf(f, "# MAPR = %d\n", header->MAPR);
	fprintf(f, "# MAPS = %d\n", header->MAPS);
	fprintf(f, "# AMIN = %f\n", header->AMIN);
	fprintf(f, "# AMAX = %e\n", header->AMAX);
	fprintf(f, "# AMEAN = %e\n", header->AMEAN);
	fprintf(f, "# ISPG = %d\n", header->ISPG);
	fprintf(f, "# NSYMBT = %d\n", header->NSYMBT);
	fprintf(f, "# LSKFLG = %d\n", header->LSKFLG);
	fprintf(f, "# SKWMAT =\n# ");
	for(i = 0; i < 9; i++) {
		fprintf(f, "%f ", header->SKWMAT[i]);
		if((i+1)%3 == 0) {
			fprintf(f, "\n# ");
		}
	}
	fprintf(f, "# SKWTRN =\n# ");
	for(i = 0; i < 3; i++) {
		fprintf(f, "%f ", header->SKWTRN[i]);
	}
	fprintf(f, "\n");
	fprintf(f, "# EXTRA =\n# ");
	for(i = 0; i < 14; i++) {
		fprintf(f, "%d ", header->EXTRA[i]);
	}
	fprintf(f, "\n");
	fprintf(f, "# MAP = %s\n", header->MAP);
	fprintf(f, "# MACHST = %d\n", header->MACHST);
	fprintf(f, "# RMS = %f\n", header->RMS);
	fprintf(f, "# NLABL = %d\n", header->NLABL);
	fprintf(f, "# LABEL_N = %s\n", header->LABEL_N);
}

void add_maps(map_data *map_1, map_data *map_2, map_data *map_out, int t_x, int t_y, int t_z)
{
	if(t_x + map_2->NX < map_1->NX) {
		map_out->NX = map_1->NX;
	} else {
		map_out->NX = t_x + map_2->NX;
	}

	if(t_y + map_2->NY < map_1->NY) {
		map_out->NY = map_1->NY;
	} else {
		map_out->NY = t_y + map_2->NY;
	}

	if(t_z + map_2->NZ < map_1->NZ) {
		map_out->NZ = map_1->NZ;
	} else {
		map_out->NZ = t_z + map_2->NZ;
	}

	printf("New map dimensions: (%d %d %d)\n", map_out->NX, map_out->NY, map_out->NZ);

	map_out->num_voxels = map_out->NX * map_out->NY * map_out->NZ;
	map_out->data = new float[map_out->num_voxels];
	for(int i = 0; i < map_out->num_voxels; i++) {
		map_out->data[i] = 0;
	}

	// add map_1 to map_out at (0,0,0)
	for(int z = 0; z < map_1->NZ; z++) {
		for(int y = 0; y < map_1->NY; y++) {
			for(int x = 0; x < map_1->NX; x++) {
				map_out->data[MAP_OUT_INDEX(x, y, z)] += map_1->data[MAP1_INDEX(x, y, z)];
			}
		}
	}

	// add map_2 to map_out at (t_x,t_y,t_z)
	for(int z = 0; z < map_2->NZ; z++) {
		for(int y = 0; y < map_2->NY; y++) {
			for(int x = 0; x < map_2->NX; x++) {
				map_out->data[MAP_OUT_INDEX(x + t_x, y + t_y, z + t_z)] += map_2->data[MAP1_INDEX(x, y, z)];
			}
		}
	}

	map_out->voxel_sizeX = map_1->voxel_sizeX;
	map_out->voxel_sizeY = map_1->voxel_sizeY;
	map_out->voxel_sizeZ = map_1->voxel_sizeZ;
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

int main(int argc, char **argv)
{
	char map_1_fname[MAX_FNAME_LENGTH], map_2_fname[MAX_FNAME_LENGTH], out_fname[MAX_FNAME_LENGTH];

	if(argc != 7) {
		printf("Usage: ./emdb_map_to_ffea INPUT_MAP_1_FNAME INPUT_MAP_2_FNAME OUTPUT_MAP_FNAME TRANS_X TRANS_Y TRANS_Z\n");
		return 1;
	}

	sprintf(map_1_fname, "%s", argv[1]);
	sprintf(map_2_fname, "%s", argv[2]);
	sprintf(out_fname, "%s", argv[3]);

	int t_x = atoi(argv[4]);
	int t_y = atoi(argv[5]);
	int t_z = atoi(argv[6]);

	// extract header and voxel data from map_1
	printf("Reading in map 1...");
	header_info header_1;
	map_data map_1;
	if(extract_data(map_1_fname, &header_1, &map_1) != 0) {
		printf("Error when trying to extract data from map file 1.\n");
		return 1;
	}
	print_header_to_file(stdout, &header_1);

	// extract header and voxel data from map_2
	printf("Reading in map 2...");
	header_info header_2;
	map_data map_2;
	if(extract_data(map_2_fname, &header_2, &map_2) != 0) {
		printf("Error when trying to extract data from map file 2.\n");
		return 1;
	}
	print_header_to_file(stdout, &header_2);

	// add map_2 to map_1 translated by (t_x, t_y, t_z)
	map_data map_result;
	add_maps(&map_1, &map_2, &map_result, t_x, t_y, t_z);

	// modify the header for map 1 and use that as map_result's header
	header_1.NC = map_result.NX;
	header_1.NR = map_result.NY;
	header_1.NS = map_result.NZ;
	header_1.NX = map_result.NX;
	header_1.NY = map_result.NY;
	header_1.NZ = map_result.NZ;
	header_1.X_LENGTH = map_result.NX * map_result.voxel_sizeX;
	header_1.Y_LENGTH = map_result.NY * map_result.voxel_sizeY;
	header_1.Z_LENGTH = map_result.NZ * map_result.voxel_sizeZ;

	// write new map to file
	write_map_to_file(&map_result, &header_1, out_fname);

	return 0;
}
