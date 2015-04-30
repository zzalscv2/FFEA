/*
 *      Blob.h
 */

#ifndef BLOB_H_INCLUDED
#define BLOB_H_INCLUDED

#include "SparsityPattern.h"
#include "ConnectivityTypes.h"

#include <stdio.h>
#include <math.h>

#include <omp.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "ConjugateGradientSolver.h"
#include "MassLumpedSolver.h"
#include "NoMassCGSolver.h"
#include "Face.h"
#include "CG_solver.h"
#include "BEM_Poisson_Boltzmann.h"
#include "LJ_matrix.h"

/*
 * The "Blob" class
 */
class Blob {
public:
    /*
     * Blob constructor:
     * Initialises all variables and pointers to 0 (or NULL). Does not perform any memory allocation. Actual Blob initialisation
     * is carried out by the init() method.
     */
    Blob();

    /*
     * Blob destructor:
     * Closes all open files, deallocates memory held by arrays, solvers etc. and sets everything to zero (NULL)
     */
    ~Blob();

    /*
     * Allocates the node, element and embedded charge arrays, initialised with the data
     * read from the given node, element and embedded charge files.
     * 'linear_solver' sets which type of linear solver is to be used: 0 for direct (forward/backward
     * substitution) and 1 for iterative (preconditioned gonjugate gradient).
     * Also takes the simulation parameters and the array of RNGs (for multiprocessor runs).
     */
    int init(const int blob_index, const int conformation_index, const char *node_filename, const char *topology_filename, const char *surface_filename, const char *material_params_filename,
            const char *stokes_filename, const char *vdw_filename, const char *pin_filename, scalar scale, int linear_solver,
            int blob_state, SimulationParams *params, LJ_matrix *lj_matrix, MTRand rng[], int num_threads);


    /*
     * Solves the EOM on the finite element mesh, updating the node positions and velocities by one time step
     */
    int update();

    /*
     * Calculates the centroid of this Blob, then brings the Blob to the origin, 
     * rotates all nodes in the Blob, and brings back the Blob to the initial position.
     * aso that the new centroid position is at the given (x,y,z) position
     */
    void rotate(float r11, float r12, float r13, float r21, float r22, float r23, float r31, float r32, float r33);

    /*
     * Calculates the centroid of this Blob, then translates all nodes in the Blob
     * so that the new centroid position is at the given (x,y,z) position
     */
    void position(scalar x, scalar y, scalar z);


    /*
     * Translate the Blob by the given vector
     */
    void move(scalar dx, scalar dy, scalar dz);

    /*
     * Calculates and returns the centre of mass of this Blob
     */
    void get_CoM(vector3 *com);

    /*
     * Calculates and returns the centroid of this Blob
     */
    void get_centroid(vector3 *com);

    void set_rmsd_pos_0();

    /*
     * Writes a new node file for the case of an initially translated STATIC blob, so viewer doesn't read straight from node file
     */
    int create_viewer_node_file(const char *node_filename, scalar scale);

    /*
     * Dumps all the node positions (in order) from the node array to the given file stream.
     */
    void write_nodes_to_file(FILE *trajectory_out);

    /*
     * Reads the node positions from the given trajectory file stream.
     * This is useful for restarting simulations from trajectory files.
     */
    int read_nodes_from_file(FILE *trajectory_out);

    /*
     * Takes measurements of system properties: KE, PE, Centre of Mass and Angular momentum, outputing
     * the results to the measurement file.
     */
    void make_measurements(FILE *measurement_out, int step, vector3 *system_CoM);

    void make_stress_measurements(FILE *stress_out, int blob_number);

    void calculate_vdw_bb_interaction_with_another_blob(FILE *vdw_measurement_out, int other_blob_index);

    /*
     * Get the centroid of all faces on the blob surface
     */
    void calc_centroids_and_normals_of_all_faces();

    /*
     *
     */
    int get_num_faces();

    /*
     * Return pointer to the ith Face of this Blob's surface
     */
    Face *get_face(int i);

    scalar get_vdw_area();

    //		/*
    //		 *
    //		 */
    //		int get_num_surface_elements();


    /*
     * Solves the poisson equation inside this Blob for a given fixed surface potential, and calculates the surface flux out of the protein
     */
    int solve_poisson(scalar *phi_gamma_IN, scalar *J_Gamma_OUT);

    /*
     * Set all forces on the blob to zero
     */
    void zero_force();

    void set_forces_to_zero();

    vector3 get_node(int index);

    void add_force_to_node(vector3 f, int index);

    void zero_vdw_bb_measurement_data();

    void zero_vdw_xz_measurement_data();

    /*
     * Set all nodes on the Blob to the given velocity vector
     */
    void velocity_all(scalar vel_x, scalar vel_y, scalar vel_z);

    /*
     * Constructs the Poisson matrix for this Blob.
     */
    void build_poisson_matrices();

    //		int elements_are_connected(int e1, int e2);

    /*
     * Returns the total mass of this Blob.
     */
    scalar get_mass();

    /*
     * Applies the WALL_TYPE_HARD boundary conditions simply by finding all nodes that have "passed through" the wall,
     * then zeroing any component of those nodes' velocities that points into the wall. This allows the Blob to "wobble"
     * its way back out of the wall (and effectively prevents any penetration larger than dt * largest velocity,
     * generally very small for sensible dt).
     */
    void enforce_box_boundaries(vector3 *box_dim);

    /*
     * Set the interaction flag for all faces of this Blob back to false (i.e. not interacting)
     */
    void reset_all_faces();

    int get_num_nodes();

    /* Blob and conformation indices */
    int blob_index;
    int conformation_index;

private:

    /* Total number of nodes in Blob */
    int num_nodes;

    /* Total number of elements in Blob */
    int num_elements;

    /* Total number of surface elements in Blob */
    int num_surface_elements;

    /* Total number of interior elements in Blob */
    int num_interior_elements;

    /* Total number of faces on Blob surface */
    int num_surface_faces;

    /* Total number of nodes on Blob surface */
    int num_surface_nodes;

    /* Total number of nodes on Blob interior */
    int num_interior_nodes;

    /* Number of 'pinned' nodes (nodes which are not able to move, removing degrees of freedom from system) */
    int num_pinned_nodes;

    /* Whether this Blob is DYNAMIC (movable; dynamics simulated) or STATIC (fixed; no simulation of dynamics; Blob is a perfectly solid object fixed in space)*/
    int blob_state;

    /* Total mass of Blob */
    scalar mass;

    /* Total vdw energy between blobs */
    scalar vdw_bb_energy;

    /* Array of nodes */
    mesh_node *node;

    /* Array of elements */
    tetra_element_linear *elem;

    /* Array of surface faces */
    Face *surface;

    /* List of fixed ('pinned') nodes */
    int *pinned_nodes_list;

    /* A pointer to a class containing simulation parameters, such as the time step, dt */
    SimulationParams *params;

    /* pointer to the vdw forcefield parameters (for some energy calcs) */
    LJ_matrix *lj_matrix;

    /* A pointer to whatever Solver is being used for this Blob (eg SparseSubstitutionSolver
     * or ConjugateGradientSolver). The Solver solves the equation Mx = f where M is the
     * mass matrix of this Blob and f is the force vector.
     */
    Solver *solver;

    /* Remember what type of solver we are using */
    int linear_solver;

    /* The Blob force vector (an array of the force on every node) */
    vector3 *force;

    /* The array of random number generators (needed for parallel runs) */
    MTRand *rng;

    CG_solver *poisson_solver;
    SparseMatrixFixedPattern *poisson_surface_matrix;
    SparseMatrixFixedPattern *poisson_interior_matrix;
    scalar *phi_Omega;
    scalar *phi_Gamma;
    scalar *q;
    scalar *nodal_q;
    scalar *poisson_rhs;

    int *num_contributing_faces;

    connectivity_entry *element_connectivity_table;

    /*
     *
     */
    SparseMatrixFixedPattern *M;


    /*
     * Opens and reads the given 'ffea node file', extracting all the nodes for this Blob.
     * Records how many of these are surface nodes and how many are interior nodes.
     */
    int load_nodes(const char *node_filename, scalar scale);

    /*
     * Opens and reads the given 'ffea topology file', extracting all the elements for this Blob.
     * Records how many of these are surface elements and how many are interior elements.
     */
    int load_topology(const char *topology_filename);

    /*
     * Opens and reads the given 'ffea surface file', extracting all the faces for this Blob.
     */
    int load_surface(const char *surface_filename, SimulationParams* params);


    /*
     * Opens and reads the given 'ffea surface file', extracting all the faces for this Blob, ignoring topology.
     */
    int load_surface_no_topology(const char *surface_filename, SimulationParams *params);

    /*
     * Opens and reads the given 'ffea material params file', extracting the material parameters for each element.
     */
    int load_material_params(const char *material_params_filename);

    /*
     * Opens and reads the given 'ffea stokes params file', extracting the stokes radii for each node in the Blob.
     */
    int load_stokes_params(const char *stokes_filename, scalar scale);


    /*
     * Opens and reads the given 'ffea vdw file', extracting all the van der waals species for each face of this Blob.
     */
    int load_vdw(const char *vdw_filename, int num_vdw_face_types);

    /*
     * Opens and reads the given 'ffea pinned nodes file', extracting all the faces for this Blob.
     */
    int load_pinned_nodes(const char *pin_filename);


    /*
     * Calculate some quantities such as the rest jacobian, rest volume etc.
     */
    void calc_rest_state_info();

    /*
     *
     */
    int aggregate_forces_and_solve();

    /*
     *
     */
    void euler_integrate();

    /*
     *
     */
    int calculate_node_element_connectivity();

    void build_mass_matrix();
};

#endif