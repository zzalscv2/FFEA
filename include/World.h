#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <typeinfo>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include "MersenneTwister.h"
#include "NearestNeighbourLinkedListCube.h"
#include "BEM_Poisson_Boltzmann.h"
#include "BiCGSTAB_solver.h"
#include "FFEA_return_codes.h"
#include "FFEA_input_reader.h"
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
#include "PreComp_solver.h"
#include "LJ_matrix.h"
#include "BindingSite.h"
#include "Spring.h"
#include "SparseMatrixFixedPattern.h"
#include "KineticState.h"
#include "KineticBindingSite.h"

#include "dimensions.h"
using namespace std;

class World {
public:
    World();

    ~World();

    /* */
    int init(string FFEA_script_filename, int frames_to_delete, int mode);

    /* */
    int get_smallest_time_constants();

    /* */
    int enm(set<int> blob_indices, int num_modes);

    /* */
    int dmm(set<int> blob_indices, int num_modes);

    /* */
    int run();

    /* */
    int read_and_build_system(vector<string> script_vector);

    /* */
    int load_kinetic_maps(vector<string> map_fnames, vector<int> map_from, vector<int> map_to, int blob_index);

    /* */
    int build_kinetic_identity_maps();

    /* */
    int load_kinetic_states(string states_fname, int blob_index);

    /* */
    int load_kinetic_rates(string rates_fname, int blob_index);

    /* */
    void get_system_CoM(vector3 *system_CoM);

    /* */
    void get_system_centroid(vector3 *centroid);

    /* */
    void get_system_dimensions(vector3 *dimenstion_vector);

    /* */
    int enm(int *blob_index, int num_modes);

    /* */
    int get_num_blobs();

private:

    /** @brief How many Blobs populate this world */
    int num_blobs;

    /** @brief How many conformations does each blob have? */
    int *num_conformations;

    /** @brief 2-D Array of Blob objects (blob i, conformation j) */
    Blob **blob_array;

    /** @brief Which conformation is active in each blob */
    Blob **active_blob_array;
    int *active_conformation_index;
    int *active_state_index;

    /** @brief Maps for kinetic switching of conformations */
    SparseMatrixFixedPattern ***kinetic_map_array;
    SparseMatrixFixedPattern ****kinetic_double_map_array;

    /** @brief Kinetic State and Rate objects */
    KineticState **kinetic_state;
    scalar ***kinetic_rate;

    /** @brief An array of springs which connect nodes if necessary */
    Spring *spring_array;

    /** @brief And how many springs are there? */
    int num_springs;

    /** @brief How many kinetic binding sites are there? */
    int num_binding_sites;

    /** @brief How many threads are available for parallelisation */
    int num_threads;

    /** @brief An array of pointers to random number generators (for use in parallel) */
    MTRand *rng;

    /** @brief An array of pointers to random number generators for use in kinetics */
    MTRand kinetic_rng;

    /** @brief Parameters being used for this simulation */
    SimulationParams params;

    /** @brief
     * Data structure keeping track of which `cell' each face lies in (where the world has been discretised into a grid of cells of dimension 1.5 kappa)
     * so that the BEM matrices may be constructed quickly and sparsely.
     */
    NearestNeighbourLinkedListCube lookup;

    /** @brief * Output trajectory file */
    FILE *trajectory_out;

    /** @brief * Output measurement file */
    FILE **measurement_out;

    /*
     *
     */
    //		SurfaceElementLookup surface_element_lookup;

    /** @brief BEM solver for the exterior electrostatics */
    BEM_Poisson_Boltzmann PB_solver;

    /** @brief Number of surface faces in entire system */
    int total_num_surface_faces;

    /** @brief
     * Vector of the electrostatic potential on each surface in entire system
     */
    scalar *phi_Gamma;
    scalar *J_Gamma;
    scalar *work_vec;

    /** @brief
     * Biconjugate gradient stabilised solver for nonsymmetric matrices
     */
    BiCGSTAB_solver nonsymmetric_solver;

    /* Van der Waals solver */
    VdW_solver vdw_solver;

    /** @brief LJ parameters matrix */
    LJ_matrix lj_matrix;

    /** @brief Binding Interactions matrix */
    BindingSite_matrix binding_matrix;


    /**
      * @brief stores info within the <precomp> block at the .ffea file.
      */
    PreComp_params pc_params; 
    /**
      * @brief PreComputed potentials solver
      */
    PreComp_solver pc_solver;


    vector3 box_dim;

    long long step_initial;

    int load_springs(const char *fname);

    void activate_springs();

    void apply_springs();

    int rescale_kinetic_rates(scalar ***rates);
	
    int change_blob_state(int blob_index, int new_state_index);

    int get_next_script_tag(FILE *in, char *buf);

    void apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N);

    void do_es();

    void make_trajectory_from_eigenvector(string traj_out_fname, int blob_index, int mode_index, Eigen::VectorXd evec, double step);

    void print_evecs_to_file(string fname, Eigen::MatrixXd ev, int num_rows, int num_modes);

    void print_evals_to_file(string fname, Eigen::VectorXd ev, int num_modes);

    void write_eig_to_files(double *evals_ordered, double **evecs_ordered, int num_modes, int num_nodes);
    
    void print_trajectory_and_measurement_files(int step, double wtime);

    void print_trajectory_conformation_changes(FILE *fout, int step, int *from_index, int *to_index);

    void print_static_trajectory(int step, double wtime, int blob_index);
};

#endif
