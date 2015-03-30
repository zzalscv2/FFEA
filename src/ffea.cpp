#include <stdio.h>
#include <stdlib.h>
#include <ctime>

#include <omp.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"

#define MAX_FNAME_LENGTH 255

int main(int argc, char *argv[]) {
    printf("\n\n\n***************************************************\n\tFLUCTUATING FINITE ELEMENT ANALYSIS\n***************************************************\n\n");
    printf(" Version:\t%s [%s]\n", FFEA_VERSION, FFEA_MASCOT);
    printf("Compiled:\t%s at %s\n", __DATE__, __TIME__);
    printf("  Coding:\tRobin Richardson (pyrar@leeds.ac.uk), Ben Hanson (py09bh@leeds.ac.uk)\n");
    printf("  Theory:\tOliver Harlen, Sarah Harris, Robin Oliver, Daniel Read, Robin Richardson, Ben Hanson\n");

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    printf("Parallelisation switch: FFEA_PARALLEL_WITHIN_BLOB\n");
#endif

#ifdef FFEA_PARALLEL_PER_BLOB
    printf("Parallelisation switch: FFEA_PARALLEL_PER_BLOB\n");
#endif

    // Return variable
    int myreturn = 0;

    // Check for correct number of arguments
    if (argc < 2) {
        printf("\n\nUsage: ffea [FFEA SCRIPT FILE (.ffea)] [OPTIONS]\n\n\n");
        return FFEA_ERROR;
    }

    // Read arguments
    //int arg = 2;

    char script_fname[MAX_FNAME_LENGTH];
    sprintf(script_fname, "%s", argv[1]);
    /*while(arg < argc - 1) {
            if(strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--mode") == 0) {
                    arg++;
                    if(strcmp(argv[arg], "s") == 0) {
                            simulation_mode = SIMULATION_MODE;
                    } else if (strcmp(argv[arg], "r") == 0) {
                            simulation_mode = RELAXATION_MODE;
                    } else {
                            FFEA_error_text();
                            printf("Accepted arguments for '-m/--mode':\n\n's'\n'r'\n\n");
                            printf("Errors during argument initialisation mean simulation cannot commence.\n");
                            return FFEA_ERROR;
                    }
            } else {
                    FFEA_error_text();
                    printf("Accepted flags are:\n\n'-m/--mode'\n\n");
                    printf("Errors during argument initialisation mean simulation cannot commence.\n");
                    return FFEA_ERROR;
            }
    }*/

    // The system of all proteins, electrostatics and water
    World *world;

    // Allocate the world
    world = new World();

    // Initialise the world, loading all blobs, parameters, electrostatics, etc.
    printf("Initialising the world:\n");
    if (world->init(script_fname) == FFEA_ERROR) {
        FFEA_error_text();
        printf("Errors during initialisation mean World cannot be constructed properly.\n");
    } else {
        // Run the world for the specified number of time steps
        printf("Running world...\n");
        if (world->run() == FFEA_ERROR) {
            FFEA_error_text();
            printf("Error occurred when running World\n");
            myreturn = FFEA_ERROR;
        } else {
            printf("...done\n");
            myreturn = FFEA_OK;
        }
    }

    // Delete the world (oh no!)
    printf("Deleting world...");
    delete world;
    printf("...done. World has been sucessfully destroyed.\n");

    return myreturn;
}
