#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <omp.h>

#include "MersenneTwister.h"
#include "NearestNeighbourLinkedListCube.h"
#include "BEM_Poisson_Boltzmann.h"
#include "BiCGSTAB_solver.h"
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
#include "VdW_solver.h"
#include "LJ_matrix.h"
#include "Spring.h"

class World {
public:
    World();

    ~World();

    /* */
    int init(const char *FFEA_script_filename);

    /* */
    int run();

    /* */
    int read_and_build_system(FILE *in);

    /* */
    void get_system_CoM(vector3 *system_CoM);

private:

    /* How many Blobs populate this world */
    int num_blobs;

    /* How many conformations does each blob have? */
    int *num_conformations;

    /* 2-D Array of Blob objects (blob i, conformation j) */
    Blob **blob_array;

    /* Which conformation is active in each blob */
    Blob **active_blob_array;
    int *active_conformation_index;

    /* An array of springs which connect nodes if necessary */
    Spring *spring_array;

    /* And how many springs are there? */
    int num_springs;

    /* How many threads are available for parallelisation */
    int num_threads;

    /* An array of pointers to random number generators (for use in parallel) */
    MTRand *rng;

    /* Parameters being used for this simulation */
    SimulationParams params;

    /*
     * Data structure keeping track of which `cell' each face lies in (where the world has been discretised into a grid of cells of dimension 1.5 kappa)
     * so that the BEM matrices may be constructed quickly and sparsely.
     */
    NearestNeighbourLinkedListCube lookup;

    /* * Output trajectory file */
    FILE *trajectory_out;

    /* * Output measurement file */
    FILE **measurement_out;

    /* * Output stress measurement file */
    FILE *stress_out;
    /*
     *
     */
    //		SurfaceElementLookup surface_element_lookup;

    /* * BEM solver for the exterior electrostatics */
    BEM_Poisson_Boltzmann PB_solver;

    /* * Number of surface faces in entire system */
    int total_num_surface_faces;

    /*
     * Vector of the electrostatic potential on each surface in entire system
     */
    scalar *phi_Gamma;
    scalar *J_Gamma;
    scalar *work_vec;

    /*
     * Biconjugate gradient stabilised solver for nonsymmetric matrices
     */
    BiCGSTAB_solver nonsymmetric_solver;

    /* * Van der Waals solver */
    VdW_solver vdw_solver;

    /* * LJ parameters matrix */
    LJ_matrix lj_matrix;

    vector3 box_dim;

    long long step_initial;

    int load_springs(char *fname);

    void activate_springs();

    void apply_springs();

    int get_next_script_tag(FILE *in, char *buf);

    void apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N);

    void do_es();

    void print_trajectory_and_measurement_files(int step, double wtime);

    void print_static_trajectory(int step, double wtime, int blob_index);
};

#endif
