#ifndef FFEA_RETURN_CODES_H_INCLUDED
#define FFEA_RETURN_CODES_H_INCLUDED

#define FFEA_VERSION "2.0"
#define FFEA_MASCOT "Mega Walrus"

#define FFEA_DIRECT_SOLVER		0
#define FFEA_ITERATIVE_SOLVER		1
#define FFEA_MASSLUMPED_SOLVER		2
#define FFEA_NOMASS_CG_SOLVER		3

#define FFEA_OK		0
#define FFEA_ERROR	-1
#define FFEA_CAUTION 	1

#define FFEA_FILE_ERROR_MESSG(F) {FFEA_error_text(); printf("Error opening file: %s\n", F); perror(NULL); return FFEA_ERROR;}
#define FFEA_ERROR_MESSG(...) {FFEA_error_text(); printf(__VA_ARGS__); return FFEA_ERROR;}

//#define FFEA_PARALLEL_WITHIN_BLOB
//#define FFEA_PARALLEL_PER_BLOB

#define FFEA_BLOB_IS_STATIC	0
#define FFEA_BLOB_IS_DYNAMIC	1
#define FFEA_BLOB_IS_FROZEN	2

#define FFEA_CONFORMATION_CHANGE	0
#define FFEA_BINDING_EVENT		1
#define FFEA_UNBINDING_EVENT		2
#define FFEA_IDENTITY_EVENT		3

#include <stdio.h>

/* Prints "ERROR: " to stderr in red text and "CAUTION: " in yellow */
void FFEA_error_text();
void FFEA_caution_text();
#endif
