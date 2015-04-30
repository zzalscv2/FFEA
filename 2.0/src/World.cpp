#include "World.h"

World::World() {
    // Initialise everything to zero
    blob_array = NULL;
    spring_array = NULL;
    num_blobs = 0;
    num_conformations = NULL;
    num_springs = 0;
    num_threads = 0;
    rng = NULL;
    phi_Gamma = NULL;
    total_num_surface_faces = 0;
    box_dim.x = 0;
    box_dim.y = 0;
    box_dim.z = 0;
    step_initial = 0;
    stress_out = NULL;
    trajectory_out = NULL;
    measurement_out = NULL;
}

World::~World() {
    delete[] rng;
    rng = NULL;
    num_threads = 0;

    delete[] blob_array;
    blob_array = NULL;
    num_blobs = 0;
    delete[] num_conformations;
    num_conformations = NULL;

    delete[] spring_array;
    spring_array = NULL;
    num_springs = 0;

    delete[] phi_Gamma;
    phi_Gamma = NULL;

    total_num_surface_faces = 0;

    box_dim.x = 0;
    box_dim.y = 0;
    box_dim.z = 0;
    step_initial = 0;

    stress_out = NULL;
    trajectory_out = NULL;

    delete[] measurement_out;
    measurement_out = NULL;
}

/* */
int World::init(const char *FFEA_script_filename) {
    int i, j;

    // Open the ffea script file
    FILE *script_file;
    if ((script_file = fopen(FFEA_script_filename, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(FFEA_script_filename);
    }

    const int buffer_size = 101;
    char buf[buffer_size];

    printf("Parsing script file '%s'...\n", FFEA_script_filename);

    // Read first tag
    if (get_next_script_tag(script_file, buf) == FFEA_ERROR) {
        FFEA_error_text();
        printf("Error when reading first tag in %s\n", FFEA_script_filename);
        fclose(script_file);
        return FFEA_ERROR;
    }

    // 1st block should be <param> block
    if (strcmp(buf, "param") != 0) {
        FFEA_error_text();
        printf("In script file, '%s': First block in script file should be <param> block, not '<%s>'\n", FFEA_script_filename, buf);
        fclose(script_file);
        return FFEA_ERROR;
    }

    // Parse each tag in the <param> block
    printf("Entering <param> block...\n");
    for (;;) {
        if (get_next_script_tag(script_file, buf) == FFEA_ERROR) {
            FFEA_error_text();
            printf("\tError when reading tag in <param> block\n");
            fclose(script_file);
            return FFEA_ERROR;
        }
        if (feof(script_file)) {
            FFEA_error_text();
            printf("\tReached end of file before end of <param> block\n");
            fclose(script_file);
            return FFEA_ERROR;
        }

        if (strcmp(buf, "/param") == 0) {
            printf("Exiting <param> block\n");
            break;
        }

        if (params.parse_param_assignment(buf) == FFEA_ERROR) {
            FFEA_error_text();
            printf("Parameter assignment failed.\n");
            fclose(script_file);
            return FFEA_ERROR;
        }
    }

    // Check if the simulation parameters are reasonable
    if (params.validate() == FFEA_ERROR) {
        FFEA_error_text();
        printf("There are some required parameters that are either unset or have bad values.\n");
        fclose(script_file);
        return FFEA_ERROR;
    }

    // Add data to world for quicker coding!
    this->num_blobs = params.num_blobs;
    this->num_conformations = new int[num_blobs];

    for (int i = 0; i < num_blobs; ++i) {
        this->num_conformations[i] = params.num_conformations[i];
    }

    // Load the vdw forcefield params matrix
    if (lj_matrix.init(params.vdw_params_fname) == FFEA_ERROR) {
        FFEA_ERROR_MESSG("Error when reading from vdw forcefeild params file.\n")
    }

    // detect how many threads we have for openmp
    int tid;
#pragma omp parallel default(none) private(tid)
    {
        tid = omp_get_thread_num();
        if (tid == 0) {
            num_threads = omp_get_num_threads();
            printf("Number of threads detected: %d\n", num_threads);
        }
    }

    // We need one rng for each thread, to avoid concurrency problems,
    // so generate an array of instances of mersenne-twister rngs.
    rng = new MTRand[num_threads];

    // Seed each rng differently
    for (i = 0; i < num_threads; i++)
        rng[i].seed(params.rng_seed + i);


    // Read next main block tag
    if (get_next_script_tag(script_file, buf) == FFEA_ERROR) {
        FFEA_error_text();
        printf("Error when reading/looking for <system> tag in %s\n", FFEA_script_filename);
        fclose(script_file);
        return FFEA_ERROR;
    }

    // 2nd block should be <system> block
    if (strcmp(buf, "system") != 0) {
        FFEA_error_text();
        printf("In script file, '%s': Second block in script file should be <system> block, not '<%s>'\n", FFEA_script_filename, buf);
        fclose(script_file);
        return FFEA_ERROR;
    }

    // Parse system block
    printf("Entering <system> block...\n");
    if (read_and_build_system(script_file) == FFEA_ERROR) {
        FFEA_error_text();
        printf("\tError when reading <system> block. Aborting.\n");
        fclose(script_file);
        return FFEA_ERROR;
    }
    fclose(script_file);

    // Create measurement files
    measurement_out = new FILE *[params.num_blobs + 1];

    // If not restarting a previous simulation, create new trajectory and measurement files
    if (params.restart == 0) {
        // Open the trajectory output file for writing
        if ((trajectory_out = fopen(params.trajectory_out_fname, "w")) == NULL) {
            FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
        }

        // Open the measurement output file for writing
        for (i = 0; i < params.num_blobs + 1; ++i) {
            if ((measurement_out[i] = fopen(params.measurement_out_fname[i], "w")) == NULL) {
                FFEA_FILE_ERROR_MESSG(params.measurement_out_fname[i])
            }
        }

        // Open stress file for writing, maybe
        if (params.stress_out_fname_set == 1) {
            if (((stress_out = fopen(params.stress_out_fname, "w")) == NULL)) {
                FFEA_FILE_ERROR_MESSG(params.stress_out_fname)
            }
        }

        // Print initial info stuff
        fprintf(trajectory_out, "FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\n", params.num_blobs);
        for (i = 0; i < params.num_blobs; ++i) {
            fprintf(trajectory_out, "Blob %d Nodes %d\t", i, active_blob_array[i]->get_num_nodes());
        }
        fprintf(trajectory_out, "\n\n");

        // First line in trajectory data should be an asterisk (used to delimit different steps for easy seek-search in restart code)
        fprintf(trajectory_out, "*\n");

        // First line in measurements file should be a header explaining what quantities are in each column
        for (i = 0; i < params.num_blobs; ++i) {
            fprintf(measurement_out[i], "# step | KE | PE | CoM x | CoM y | CoM z | L_x | L_y | L_z | rmsd | vdw_area_%d_surface | vdw_force_%d_surface | vdw_energy_%d_surface\n", i, i, i);
            fflush(measurement_out[i]);
        }
        fprintf(measurement_out[params.num_blobs], "# step ");
        for (i = 0; i < params.num_blobs; ++i) {
            for (j = i + 1; j < params.num_blobs; ++j) {
                fprintf(measurement_out[params.num_blobs], "| vdw_area_%d_%d | vdw_force_%d_%d | vdw_energy_%d_%d ", i, j, i, j, i, j);
            }
        }
        fprintf(measurement_out[params.num_blobs], "\n");
        fflush(measurement_out[params.num_blobs]);

    } else {
        // Otherwise, seek backwards from the end of the trajectory file looking for '*' character (delimitter for snapshots)
        printf("Restarting from trajectory file %s\n", params.trajectory_out_fname);
        if ((trajectory_out = fopen(params.trajectory_out_fname, "r")) == NULL) {
            FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
        }

        printf("Reverse searching for 2 asterisks (denoting a completely written snapshot)...\n");
        if (fseek(trajectory_out, 0, SEEK_END) != 0) {
            FFEA_ERROR_MESSG("Could not seek to end of file\n")
        }

        // Variable to store position of last asterisk in trajectory file (initialise it at end of file)
        off_t last_asterisk_pos = ftello(trajectory_out);

        int num_asterisks = 0;
        while (num_asterisks != 2) {
            if (fseek(trajectory_out, -2, SEEK_CUR) != 0) {
                perror(NULL);
                FFEA_ERROR_MESSG("Balls.\n")
            }
            char c = fgetc(trajectory_out);
            if (c == '*') {
                num_asterisks++;
                printf("Found %d\n", num_asterisks);

                // get the position in the file of this last asterisk
                if (num_asterisks == 1) {
                    last_asterisk_pos = ftello(trajectory_out);
                }
            }
        }

        char c;
        if ((c = fgetc(trajectory_out)) != '\n') {
            ungetc(c, trajectory_out);
        }

        printf("Loading Blob position and velocity data from last completely written snapshot \n");
        int blob_id, conformation_id;
        long long rstep;
        for (int b = 0; b < params.num_blobs; b++) {
            if (fscanf(trajectory_out, "Blob %d, Conformation %d, step %lld\n", &blob_id, &conformation_id, &rstep) != 3) {
                FFEA_ERROR_MESSG("Error reading header info for Blob %d\n", b)
            }
            if (blob_id != b) {
                FFEA_ERROR_MESSG("Mismatch in trajectory file - found blob id = %d, was expecting blob id = %d\n", blob_id, b)
            }
            printf("Loading node position, velocity and potential from restart trajectory file, for blob %d, step %lld\n", blob_id, rstep);
            if (active_blob_array[b]->read_nodes_from_file(trajectory_out) == FFEA_ERROR) {
                FFEA_ERROR_MESSG("Error restarting blob %d\n", blob_id)
            }
        }
        step_initial = rstep;
        printf("...done. Simulation will commence from step %lld\n", step_initial);
        fclose(trajectory_out);

        // Truncate the trajectory file up to the point of the last asterisk (thereby erasing any half-written time steps that may occur after it)
        printf("Truncating the trajectory file to the last asterisk (to remove any half-written time steps)\n");
        if (truncate(params.trajectory_out_fname, last_asterisk_pos) != 0) {
            FFEA_ERROR_MESSG("Error when trying to truncate trajectory file %s\n", params.trajectory_out_fname)
        }

        // Open trajectory and measurment files for appending
        printf("Opening trajectory and measurement files for appending.\n");
        if ((trajectory_out = fopen(params.trajectory_out_fname, "a")) == NULL) {
            FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
        }
        for (i = 0; i < params.num_blobs + 1; ++i) {
            if ((measurement_out[i] = fopen(params.measurement_out_fname[i], "a")) == NULL) {
                FFEA_FILE_ERROR_MESSG(params.measurement_out_fname[i])
            }
        }

        // Open stress file for appending
        if (params.stress_out_fname_set == 1) {
            if (((stress_out = fopen(params.stress_out_fname, "a")) == NULL)) {
                FFEA_FILE_ERROR_MESSG(params.stress_out_fname)
            }
        }

        // Append a newline to the end of this truncated trajectory file (to replace the one that may or may not have been there)
        fprintf(trajectory_out, "\n");

        for (i = 0; i < params.num_blobs + 1; ++i) {
            fprintf(measurement_out[i], "#==RESTART==\n");
        }

    }



    // Initialise the Van der Waals solver
    box_dim.x = params.es_h * (1.0 / params.kappa) * params.es_N_x;
    box_dim.y = params.es_h * (1.0 / params.kappa) * params.es_N_y;
    box_dim.z = params.es_h * (1.0 / params.kappa) * params.es_N_z;
    vdw_solver.init(&lookup, &box_dim, &lj_matrix);

    // Calculate the total number of vdw interacting faces in the entire system
    total_num_surface_faces = 0;
    for (i = 0; i < params.num_blobs; i++) {
        total_num_surface_faces += active_blob_array[i]->get_num_faces();
    }
    printf("Total number of surface faces in system: %d\n", total_num_surface_faces);

    if (params.es_N_x > 0 && params.es_N_y > 0 && params.es_N_z > 0) {

        // Allocate memory for an NxNxN grid, with a 'pool' for the required number of surface faces
        printf("Allocating memory for nearest neighbour lookup grid...\n");
        if (lookup.alloc(params.es_N_x, params.es_N_y, params.es_N_z, total_num_surface_faces) == FFEA_ERROR) {
            FFEA_error_text();
            printf("When allocating memory for nearest neighbour lookup grid\n");
            return FFEA_ERROR;
        }
        printf("...done\n");

        printf("Box has volume %e cubic angstroms\n", (box_dim.x * box_dim.y * box_dim.z) * 1e30);

        // Add all the faces from each Blob to the lookup pool
        printf("Adding all faces to nearest neighbour grid lookup pool\n");
        for (i = 0; i < params.num_blobs; i++) {
            int num_faces_added = 0;
            for (j = 0; j < active_blob_array[i]->get_num_faces(); j++) {
                Face *b_face = active_blob_array[i]->get_face(j);
                if (b_face != NULL) {
                    if (lookup.add_to_pool(b_face) == FFEA_ERROR) {
                        FFEA_error_text();
                        printf("When attempting to add a face to the lookup pool\n");
                        return FFEA_ERROR;
                    }
                    num_faces_added++;
                }
            }
            printf("%d 'VdW active' faces, from blob %d, added to lookup grid.\n", num_faces_added, i);
        }

        // Initialise the BEM PBE solver
        if (params.calc_es == 1) {
            printf("Initialising Boundary Element Poisson Boltzmann Solver...\n");
            PB_solver.init(&lookup);
            PB_solver.set_kappa(params.kappa);

            // Initialise the nonsymmetric matrix solver
            nonsymmetric_solver.init(total_num_surface_faces, params.epsilon2, params.max_iterations_cg);

            // Allocate memory for surface potential vector
            phi_Gamma = new scalar[total_num_surface_faces];
            for (i = 0; i < total_num_surface_faces; i++)
                phi_Gamma[i] = 0;

            work_vec = new scalar[total_num_surface_faces];
            for (i = 0; i < total_num_surface_faces; i++)
                work_vec[i] = 0;

            J_Gamma = new scalar[total_num_surface_faces];
            for (i = 0; i < total_num_surface_faces; i++)
                J_Gamma[i] = 0;
        }
    }

    if (params.restart == 0) {
        // Carry out measurements on the system before carrying out any updates (if this is step 0)
        print_trajectory_and_measurement_files(0, 0);
    }

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    printf("Now ready to run with 'within-blob parallelisation' (FFEA_PARALLEL_WITHIN_BLOB) on %d threads.\n", num_threads);
#endif

#ifdef FFEA_PARALLEL_PER_BLOB
    printf("Now ready to run with 'per-blob parallelisation' (FFEA_PARALLEL_PER_BLOB) on %d threads.\n", num_threads);
#endif

    return FFEA_OK;
}

/* */
int World::run() {
    // Update entire World for num_steps time steps
    int es_count = 1;
    double wtime = omp_get_wtime();
    for (long long step = step_initial; step < params.num_steps; step++) {

        // Zero the force across all blobs
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(runtime)
#endif
        for (int i = 0; i < params.num_blobs; i++) {
            active_blob_array[i]->zero_force();
            if (params.calc_vdw == 1) {
                active_blob_array[i]->zero_vdw_bb_measurement_data();
            }
            if (params.sticky_wall_xz == 1) {
                active_blob_array[i]->zero_vdw_xz_measurement_data();
            }
        }

#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(guided) shared(stderr)
#endif
        for (int i = 0; i < params.num_blobs; i++) {
	
            // If blob centre of mass moves outside simulation box, apply PBC to it
            vector3 com;
            active_blob_array[i]->get_centroid(&com);
            scalar dx = 0, dy = 0, dz = 0;
            int check_move = 0;

            if (com.x < 0) {
                if (params.wall_x_1 == WALL_TYPE_PBC) {
                    dx += box_dim.x;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            } else if (com.x > box_dim.x) {
                if (params.wall_x_2 == WALL_TYPE_PBC) {
                    dx -= box_dim.x;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            }
            if (com.y < 0) {
                if (params.wall_y_1 == WALL_TYPE_PBC) {
                    dy += box_dim.y;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            } else if (com.y > box_dim.y) {
                if (params.wall_y_2 == WALL_TYPE_PBC) {
                    dy -= box_dim.y;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            }
            if (com.z < 0) {
                if (params.wall_z_1 == WALL_TYPE_PBC) {
                    dz += box_dim.z;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            } else if (com.z > box_dim.z) {
                if (params.wall_z_2 == WALL_TYPE_PBC) {
                    dz -= box_dim.z;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            }
            if (check_move == 1) {
                active_blob_array[i]->move(dx, dy, dz);
            }

            // If Blob is near a hard wall, prevent it from moving further into it
            active_blob_array[i]->enforce_box_boundaries(&box_dim);
        }



        if (params.calc_es == 1 || params.calc_vdw == 1 || params.sticky_wall_xz == 1) {
            if (es_count == params.es_update) {


#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(guided)
#endif
                for (int i = 0; i < params.num_blobs; i++) {
                    active_blob_array[i]->calc_centroids_and_normals_of_all_faces();
                    active_blob_array[i]->reset_all_faces();
                }


                // Attempt to place all faces in the nearest neighbour lookup table
                if (lookup.build_nearest_neighbour_lookup(params.es_h * (1.0 / params.kappa)) == FFEA_ERROR) {
                    FFEA_error_text();
                    printf("When trying to place faces in nearest neighbour lookup table.\n");

                    // attempt to print out the final (bad) time step
                    printf("Dumping final step:\n");
                    print_trajectory_and_measurement_files(step, wtime);

                    return FFEA_ERROR;
                }

                if (params.calc_vdw == 1) {
                    vdw_solver.solve();
                }
                if (params.sticky_wall_xz == 1) {
                    vdw_solver.solve_sticky_wall(params.es_h * (1.0 / params.kappa));
                }

                if (params.calc_es == 1) {
                    do_es();
                }

                es_count = 1;
            } else
                es_count++;
        }

        // Update all Blobs in the World

        // Set node forces to zero
        for (int i = 0; i < params.num_blobs; i++) {
            active_blob_array[i]->set_forces_to_zero();
        }

        // Apply springs directly to nodes
        apply_springs();

        // Sort internal forces out
        int fatal_errors = 0;
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) shared(step, wtime) reduction(+: fatal_errors) schedule(runtime)
#endif
        for (int i = 0; i < params.num_blobs; i++) {
            if (active_blob_array[i]->update() == FFEA_ERROR) {
                FFEA_error_text();
                printf("A problem occurred when updating Blob %d on step %lld\n", i, step);
                printf("Simulation ran for %2f seconds (wall clock time) before error ocurred\n", (omp_get_wtime() - wtime));
                //return FFEA_ERROR;
                fatal_errors++;
            }
        }
        if (fatal_errors > 0) {
            FFEA_error_text();
            printf("Detected %d fatal errors in this system update. Exiting now...\n", fatal_errors);

            // attempt to print out the final (bad) time step
            printf("Dumping final step:\n");
            print_trajectory_and_measurement_files(step, wtime);

            return FFEA_ERROR;
        }

        if (step % params.check == 0) {
            print_trajectory_and_measurement_files(step + 1, wtime);
        }
    }
    printf("Time taken: %2f seconds\n", (omp_get_wtime() - wtime));

    return FFEA_OK;
}

/* */
int World::read_and_build_system(FILE *in) {
    const int max_buf_size = 255;
    char buf[max_buf_size];
    char lvalue[max_buf_size];
    char rvalue[max_buf_size];
    char *spring_filename = new char[max_buf_size];
    int linear_solver = FFEA_ITERATIVE_SOLVER;
    float com_x, com_y, com_z;
    float r11, r12, r13, r21, r22, r23, r31, r32, r33;
    float vel_x, vel_y, vel_z;

    int set_linear_solver = 0;
    int set_centroid_pos = 0;
    int set_rotation = 0;
    int set_velocity = 0;

    scalar scale = 1;

    // Create the array of Blobs and conformations
    // Default active conformation is always first specified
    printf("\tCreating blob array\n");
    blob_array = new Blob*[params.num_blobs];
    active_blob_array = new Blob*[params.num_blobs];
    active_conformation_index = new int[params.num_blobs];

    for (int i = 0; i < params.num_blobs; ++i) {
        blob_array[i] = new Blob[num_conformations[i]];
        active_blob_array[i] = &blob_array[i][0];
        active_conformation_index[i] = 0;
    }

    // Construct each blob and conformation defined in the <system> block
    int i = 0, j, k = 0, rv;
    for (;;) {

        // Allocate all memory for this particular blob
        char **node_filename = new char*[num_conformations[i]];
        char **topology_filename = new char*[num_conformations[i]];
        char **surface_filename = new char*[num_conformations[i]];
        char **material_params_filename = new char*[num_conformations[i]];
        char **stokes_filename = new char*[num_conformations[i]];
        char **vdw_filename = new char*[num_conformations[i]];
        char **pin_filename = new char*[num_conformations[i]];
        int *blob_motion_state = new int[num_conformations[i]];
        char rates_filename[max_buf_size];
        char states_filename[max_buf_size];
        char map_filename[max_buf_size];

        int *set_node_filename = new int[num_conformations[i]];
        int *set_topology_filename = new int[num_conformations[i]];
        int *set_surface_filename = new int[num_conformations[i]];
        int *set_material_params_filename = new int[num_conformations[i]];
        int *set_stokes_filename = new int[num_conformations[i]];
        int *set_vdw_filename = new int[num_conformations[i]];
        int *set_pin_filename = new int[num_conformations[i]];
        int *set_blob_motion_state = new int[num_conformations[i]];
        int set_rates_filename = 0;
        int set_states_filename = 0;
        int set_map_filename = 0;

        for (j = 0; j < num_conformations[i]; ++j) {
            node_filename[j] = new char[max_buf_size];
            topology_filename[j] = new char[max_buf_size];
            surface_filename[j] = new char[max_buf_size];
            material_params_filename[j] = new char[max_buf_size];
            stokes_filename[j] = new char[max_buf_size];
            vdw_filename[j] = new char[max_buf_size];
            pin_filename[j] = new char[max_buf_size];
            blob_motion_state[j] = FFEA_BLOB_IS_DYNAMIC;

            set_node_filename[j] = 0;
            set_topology_filename[j] = 0;
            set_surface_filename[j] = 0;
            set_material_params_filename[j] = 0;
            set_stokes_filename[j] = 0;
            set_vdw_filename[j] = 0;
            set_pin_filename[j] = 0;
            set_blob_motion_state[j] = 0;

        }
        if (get_next_script_tag(in, buf) == FFEA_ERROR) {
            FFEA_error_text();
            printf("\tError reading tag in <system> block\n");
            return FFEA_ERROR;
        }

        // Check if we have prematurely reached the end of the <system> block
        if (strcmp(buf, "/system") == 0) {
            printf("\tExiting <system> block.\n");

            // Free all memory!
            delete[] node_filename;
            node_filename = NULL;
            delete[] topology_filename;
            topology_filename = NULL;
            delete[] surface_filename;
            surface_filename = NULL;
            delete[] material_params_filename;
            material_params_filename = NULL;
            delete[] stokes_filename;
            stokes_filename = NULL;
            delete[] vdw_filename;
            vdw_filename = NULL;
            delete[] pin_filename;
            pin_filename = NULL;
            delete[] blob_motion_state;
            blob_motion_state = NULL;

            delete[] set_node_filename;
            set_node_filename = NULL;
            delete[] set_topology_filename;
            set_topology_filename = NULL;
            delete[] set_surface_filename;
            set_surface_filename = NULL;
            delete[] set_material_params_filename;
            set_material_params_filename = NULL;
            delete[] set_stokes_filename;
            set_stokes_filename = NULL;
            delete[] set_vdw_filename;
            set_vdw_filename = NULL;
            delete[] set_pin_filename;
            set_pin_filename = NULL;
            delete[] set_blob_motion_state;
            set_blob_motion_state = NULL;

            break;
        }

        // If it's a blob block
        if (strcmp(buf, "blob") == 0) {

            // Read all blob initialisation info
            printf("\tEntering <blob> block %d...\n", i);

            // Get next tag (could be conformation, switching, solver, scale, centroid_pos or velocity)
            j = 0;
            for (;;) {

                // Check for errors first
                if (get_next_script_tag(in, buf) == FFEA_ERROR) {
                    FFEA_error_text();
                    printf("\t\tError reading tag in <blob> block\n");
                    return FFEA_ERROR;
                }

                if (feof(in)) {
                    FFEA_error_text();
                    printf("\t\tReached end of file before end of <blob> block\n");
                    return FFEA_ERROR;
                }

                if (strcmp(buf, "/system") == 0) {
                    FFEA_error_text();
                    printf("\t\t</system> came before </blob>\n");
                    return FFEA_ERROR;
                }

                if (strcmp(buf, "/blob") == 0) {
                    printf("\tExiting <blob> block %d\n", i);
                    break;
                }

                // Now check for possible input values
                if (strcmp(buf, "conformation") == 0) {
                    printf("\tEntering conformation block %d...\n", j);

                    // Read in node, topology, surface etc files for conformation
                    for (;;) {

                        if (get_next_script_tag(in, buf) == FFEA_ERROR) {
                            FFEA_error_text();
                            printf("\t\tError reading tag in <conformation> block\n");
                            return FFEA_ERROR;
                        }

                        if (feof(in)) {
                            FFEA_error_text();
                            printf("\t\tReached end of file before end of <conformation> block\n");
                            return FFEA_ERROR;
                        }

                        if (strcmp(buf, "/system") == 0) {
                            FFEA_error_text();
                            printf("\t\t</system> came before </conformation>\n");
                            return FFEA_ERROR;
                        }

                        if (strcmp(buf, "/blob") == 0) {
                            FFEA_error_text();
                            printf("\t\t</blob> came before </conformation>\n");
                            return FFEA_ERROR;
                        }

                        if (strcmp(buf, "/conformation") == 0) {
                            printf("\tExiting <conformation> block %d\n", j++);
                            break;
                        }

                        rv = sscanf(buf, "%100[^=]=%s", lvalue, rvalue);
                        if (rv != 2) {
                            FFEA_error_text();
                            printf("\t\tError parsing conformation parameter assignment, '%s'\n", buf);
                            return FFEA_ERROR;
                        }

                        rv = sscanf(lvalue, "%s", lvalue);
                        rv = sscanf(rvalue, "%s", rvalue);

                        if (strcmp(lvalue, "nodes") == 0) {
                            strcpy(node_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, nodes input file = '%s'\n", i, j, node_filename[j]);
                            set_node_filename[j] = 1;
                        } else if (strcmp(lvalue, "topology") == 0) {
                            strcpy(topology_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, topology input file = '%s'\n", i, j, topology_filename[j]);
                            set_topology_filename[j] = 1;
                        } else if (strcmp(lvalue, "surface") == 0) {
                            strcpy(surface_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, surface input file = '%s'\n", i, j, surface_filename[j]);
                            set_surface_filename[j] = 1;
                        } else if (strcmp(lvalue, "material") == 0) {
                            strcpy(material_params_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, material parameters file = '%s'\n", i, j, material_params_filename[j]);
                            set_material_params_filename[j] = 1;
                        } else if (strcmp(lvalue, "stokes") == 0) {
                            strcpy(stokes_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, stokes file = '%s'\n", i, j, stokes_filename[j]);
                            set_stokes_filename[j] = 1;
                        } else if (strcmp(lvalue, "vdw") == 0) {
                            strcpy(vdw_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, vdw file = '%s'\n", i, j, vdw_filename[j]);
                            set_vdw_filename[j] = 1;
                        } else if (strcmp(lvalue, "pin") == 0) {
                            strcpy(pin_filename[j], rvalue);
                            printf("\t\tSetting blob %d, conformation %d, pin file = '%s'\n", i, j, pin_filename[j]);
                            set_pin_filename[j] = 1;
                        } else if (strcmp(lvalue, "motion_state") == 0) {
                            if (strcmp(rvalue, "STATIC") == 0) {
                                printf("\t\tSetting blob %d, conformation %d, motion_state = STATIC (Fixed position, no dynamics simulated)\n", i, j);
                                blob_motion_state[j] = FFEA_BLOB_IS_STATIC;
                            } else if (strcmp(rvalue, "DYNAMIC") == 0) {
                                printf("\t\tSetting blob %d, conformation %d, motion_state = DYNAMIC (Dynamics will be simulated)\n", i, j);
                                blob_motion_state[j] = FFEA_BLOB_IS_DYNAMIC;
                            } else if (strcmp(rvalue, "FROZEN") == 0) {
                                printf("\t\tSetting blob %d, conformation %d, motion_state = FROZEN (Dynamics will not be simulated, but Blob positions still written to trajectory)\n", i, j);
                                blob_motion_state[j] = FFEA_BLOB_IS_FROZEN;
                            } else {
                                FFEA_error_text();
                                printf("\t\tFor 'motion_state' in blob %d, conformation %d, There is no state option '%s'. Please use 'STATIC', 'DYNAMIC' or 'FROZEN'.\n", i, j, rvalue);
                                return FFEA_ERROR;
                            }
                            set_blob_motion_state[j] = 1;
                        } else {
                            FFEA_error_text();
                            printf("\t\tError: In blob %d, conformation %d, '%s' is not a recognised lvalue\n", i, j, lvalue);
                            printf("\t\tRecognised lvalues are:\n");
                            printf("\t\tnodes\n\ttopology\n\tsurface\n\tmaterial\n\tstokes\n\tvdw\n\tpin\n\tmotion_state\n");
                            return FFEA_ERROR;
                        }
                    }
                } else if (strcmp(buf, "switching") == 0) {

                    // Getting conformational switching parameters map, rates and states
                    printf("\tEntering switching block...\n");
                    for (;;) {
                        if (get_next_script_tag(in, buf) == FFEA_ERROR) {
                            FFEA_error_text();
                            printf("\t\tError reading tag in <switching> block\n");
                            return FFEA_ERROR;
                        }

                        if (feof(in)) {
                            FFEA_error_text();
                            printf("\t\tReached end of file before end of <switching> block\n");
                            return FFEA_ERROR;
                        }

                        if (strcmp(buf, "/system") == 0) {
                            FFEA_error_text();
                            printf("\t\t</system> came before </switching>\n");
                            return FFEA_ERROR;
                        }
                        if (strcmp(buf, "/blob") == 0) {
                            FFEA_error_text();
                            printf("\t\t</blob> came before </switching>\n");
                            return FFEA_ERROR;
                        }
                        if (strcmp(buf, "/conformation") == 0) {
                            FFEA_error_text();
                            printf("\t\t</conformation> came before </switching>\n");
                            return FFEA_ERROR;
                        }
                        if (strcmp(buf, "/switching") == 0) {
                            printf("\tExiting <switching> block\n");
                            break;
                        }

                        rv = sscanf(buf, "%100[^=]=%s", lvalue, rvalue);
                        if (rv != 2) {
                            FFEA_error_text();
                            printf("\t\tError parsing conformation parameter assignment, '%s'\n", buf);
                            return FFEA_ERROR;
                        }

                        rv = sscanf(lvalue, "%s", lvalue);
                        rv = sscanf(rvalue, "%s", rvalue);
                        if (strcmp(lvalue, "map") == 0) {
                            strcpy(map_filename, rvalue);
                            printf("\t\tSetting blob %d, map file = '%s'\n", i, map_filename);
                            set_map_filename = 1;
                        } else if (strcmp(lvalue, "rates") == 0) {
                            strcpy(rates_filename, rvalue);
                            printf("\t\tSetting blob %d, rates file = '%s'\n", i, rates_filename);
                            set_rates_filename = 1;
                        } else if (strcmp(lvalue, "states") == 0) {
                            strcpy(states_filename, rvalue);
                            printf("\t\tSetting blob %d, states file = '%s'\n", i, states_filename);
                            set_states_filename = 1;
                        } else {
                            FFEA_error_text();
                            printf("\t\tError: In blob %d, switching block. '%s' is not a recognised lvalue\n", i, lvalue);
                            printf("\t\tRecognised lvalues are:\n");
                            printf("\t\tmap\n\trates\n\tstates\n");
                            return FFEA_ERROR;
                        }
                    }
                } else {

                    // Tag did not open a new block, so it must be a variable
                    rv = sscanf(buf, "%100[^=]=%s", lvalue, rvalue);
                    if (rv != 2) {
                        FFEA_error_text();
                        printf("\t\tError parsing conformation parameter assignment, '%s'\n", buf);
                        return FFEA_ERROR;
                    }

                    rv = sscanf(lvalue, "%s", lvalue);
                    rv = sscanf(rvalue, "%s", rvalue);
                    if (strcmp(lvalue, "solver") == 0) {
                        if (strcmp(rvalue, "CG") == 0) {
                            printf("\t\tSetting blob %d, linear solver = CG (Preconditioned Jacobi Conjugate Gradient)\n", i);
                            linear_solver = FFEA_ITERATIVE_SOLVER;
                        } else if (strcmp(rvalue, "direct") == 0) {
                            printf("\t\tSetting blob %d, linear solver = direct (Sparse Forward/Backward Substitution)\n", i);
                            linear_solver = FFEA_DIRECT_SOLVER;
                        } else if (strcmp(rvalue, "masslumped") == 0) {
                            printf("\t\tSetting blob %d, linear solver = masslumped (diagonal mass matrix)\n", i);
                            linear_solver = FFEA_MASSLUMPED_SOLVER;
                        } else if (strcmp(rvalue, "CG_nomass") == 0) {
                            printf("\t\tSetting blob %d, linear solver = CG_nomass (Preconditioned Jacobi Conjugate Gradient, no mass within system)\n", i);
                            linear_solver = FFEA_NOMASS_CG_SOLVER;
                        } else {
                            FFEA_error_text();
                            printf("\t\tThere is no solver option '%s'. Please use 'CG', 'direct' or 'CG_nomass'.\n", rvalue);
                            return FFEA_ERROR;
                        }
                        set_linear_solver = 1;
                    } else if (strcmp(lvalue, "scale") == 0) {
                        scale = atof(rvalue);
                        printf("\t\tSetting blob %d, scale factor = %e\n", i, scale);
                    } else if (strcmp(lvalue, "centroid_pos") == 0) {
                        if (sscanf(rvalue, "(%e,%e,%e)", &com_x, &com_y, &com_z) != 3) {
                            FFEA_error_text();
                            printf("\t\tCould not carry out centroid position assignment for blob %d: rvalue '%s' is badly formed. Should have form (x,y,z)\n", i, rvalue);
                            return FFEA_ERROR;
                        } else {
                            printf("\t\tSetting blob %d, initial centroid (%e, %e, %e)\n", i, com_x, com_y, com_z);
                        }
                        set_centroid_pos = 1;
                    } else if(strcmp(lvalue, "rotation") == 0)  {
                        if(sscanf(rvalue, "(%e,%e,%e,%e,%e,%e,%e,%e,%e)", &r11, &r12, &r13, &r21, &r22, &r23, &r31, &r32, &r33) != 9) {
                           FFEA_error_text();
                           printf("\t\tCould not read rotation for blob %d: rvalue '%s' is badly formed. Should have form (r11,r12,r13,r21,r22,r23,r31,r32,r33)\n", i, rvalue);
                           return FFEA_ERROR;
                        } else {
                           printf("\t\tSetting blob %d, initial rotation (%e, %e, %e, %e, %e, %e, %e, %e, %e)\n", i, r11, r12, r13, r21, r22, r23, r31, r32, r33);
                        }
                        set_rotation = 1;
                    } else if (strcmp(lvalue, "velocity") == 0) {
                        if (sscanf(rvalue, "(%e,%e,%e)", &vel_x, &vel_y, &vel_z) != 3) {
                            FFEA_error_text();
                            printf("\t\tCould not carry out velocity assignment for blob %d: rvalue '%s' is badly formed. Should have form (velx,vely,velz)\n", i, rvalue);
                            return FFEA_ERROR;
                        } else {
                            printf("\t\tSetting blob %d, initial velocity (%e, %e, %e)\n", i, vel_x, vel_y, vel_z);
                        }
                        set_velocity = 1;
                    } else if (strcmp(lvalue, "initial_conformation") == 0) {
                        if (atoi(rvalue) < 0 || atoi(rvalue) >= num_conformations[i]) {
                            FFEA_error_text();
                            printf("\t\tIn blob %d, specified 'initial_conformation' was %d. Cannot be less than 0 or greater than %d (indexing starts at 0)\n", i, atoi(rvalue), num_conformations[i]);
                            return FFEA_ERROR;
                        }
                        active_blob_array[i] = &blob_array[i][atoi(rvalue)];
                    } else {
                        FFEA_error_text();
                        printf("\t\tError: In blob %d, '%s' is not a recognised lvalue\n", i, lvalue);
                        printf("\t\tRecognised tags are:\n");
                        printf("\t\tconformation\n\tswitching\n");
                        printf("\t\tRecognised lvalues are:\n");
                        printf("\t\tscale\n\tcentroid_pos\n\tvelocity\n\tsolver\n\tinitial_conformation\n");
                        return FFEA_ERROR;
                    }
                }
            }

            // First make sure we have all the information necessary for initialisation
            if (num_conformations[i] > 1) {
                if (set_rates_filename == 0) {
                    FFEA_ERROR_MESSG("\tIn blob %d, switching block, 'rates' has not been specified\n", i);
                }
                if (set_states_filename == 0) {
                    FFEA_ERROR_MESSG("\tIn blob %d, switching block, 'states' has not been specified\n", i);
                }
                if (set_map_filename == 0) {
                    FFEA_ERROR_MESSG("\tIn blob %d, switching block, 'map' has not been specified\n", i);
                }
            }

            for (j = 0; j < num_conformations[i]; ++j) {
                if (set_blob_motion_state[j] == 0) {
                    printf("In blob %d, conformation %d: 'motion_state' option unset. Using 'DYNAMIC' by default.\n", i, j);
                    blob_motion_state[j] = FFEA_BLOB_IS_DYNAMIC;
                }

                if (set_node_filename[j] == 0) {
                    FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: 'nodes' has not been specified\n", i, j);
                }
                if (set_topology_filename[j] == 0) {
                    if (blob_motion_state[j] == FFEA_BLOB_IS_DYNAMIC) {
                        FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: 'topology' has not been specified\n", i, j);
                    }
                }
                if (set_surface_filename[j] == 0) {
                    FFEA_error_text();
                    printf("\tIn blob %d, conformation %d: 'surface' has not been specified\n", i, j);
                    return FFEA_ERROR;
                }

                if (blob_motion_state[j] == FFEA_BLOB_IS_DYNAMIC) {
                    if (set_material_params_filename[j] == 0) {
                        FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: material params file, 'material', has not been specified\n", i, j)
                    }
                    if (set_stokes_filename[j] == 0) {
                        FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: stokes radii file, 'stokes', has not been specified\n", i, j)
                    }
                    if (set_linear_solver == 0) {
                        FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: 'solver' has not been specified\n", i, j);
                    }
                }
                if (set_vdw_filename[j] == 0) {
                    FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: vdw parameters file, 'vdw', has not been specified\n", i, j)
                }
                if (set_pin_filename[j] == 0) {
                    if (blob_motion_state[j] == FFEA_BLOB_IS_DYNAMIC) {
                        FFEA_ERROR_MESSG("\tIn blob %d, conformation %d: pinned nodes file, 'pin', has not been specified\n", i, j)
                    }
                }

                printf("\tInitialising blob %d conformation %d...\n", i, j);
                if (blob_array[i][j].init(i, j, node_filename[j], topology_filename[j], surface_filename[j], material_params_filename[j], stokes_filename[j], vdw_filename[j], pin_filename[j],
                        scale, linear_solver, blob_motion_state[j], &params, &lj_matrix, rng, num_threads) == FFEA_ERROR) {
                    FFEA_error_text();
                    printf("\tError when trying to initialise Blob %d, conformation %d.\n", i, j);
                    return FFEA_ERROR;
                }

                // if centroid position is set, position the blob's centroid at that position. If vdw is set, move to center of box
                if (set_centroid_pos == 1) {
                    // Rescale first	
                    com_x *= scale;
                    com_y *= scale;
                    com_z *= scale;
                    blob_array[i][j].position(com_x, com_y, com_z);
                }

                if(set_rotation == 1) {
                    blob_array[i][j].rotate(r11,r12,r13,r21,r22,r23,r31,r32,r33);
                }

                if (set_velocity == 1)
                    blob_array[i][j].velocity_all(vel_x, vel_y, vel_z);

                // set the current node positions as pos_0 for this blob, so that all rmsd values
                // are calculated relative to this conformation centred at this point in space.
                blob_array[i][j].set_rmsd_pos_0();
            }
            i++;
        } else if (strcmp(buf, "spring") == 0) {
            for (;;) {
                if (get_next_script_tag(in, buf) == FFEA_ERROR) {
                    FFEA_error_text();
                    printf("\t\tError reading tag in <spring> block\n");
                    return FFEA_ERROR;
                }

                if (feof(in)) {
                    FFEA_error_text();
                    printf("\t\tReached end of file before end of <spring> block\n");
                    return FFEA_ERROR;
                }
                if (strcmp(buf, "/system") == 0) {
                    FFEA_error_text();
                    printf("\t\t</system> came before </spring>\n");
                    return FFEA_ERROR;
                }

                if (strcmp(buf, "/spring") == 0) {
                    printf("\t\tExiting <spring> block\n");
                    break;
                }

                // Only acceptable input is spring filename!!! 
                // By the way, I hope you're enjoying reading this source code :)
                rv = sscanf(buf, "%100[^=]=%s", lvalue, rvalue);
                if (rv != 2) {
                    FFEA_error_text();
                    printf("\t\tError parsing conformation parameter assignment, '%s'\n", buf);
                    return FFEA_ERROR;
                }

                rv = sscanf(lvalue, "%s", lvalue);
                rv = sscanf(rvalue, "%s", rvalue);
                if (strcmp(lvalue, "spring_fname") == 0) {
                    strcpy(spring_filename, rvalue);
                    printf("\t\tSetting springs input file = '%s'\n", spring_filename);
                    k++;
                } else {
                    FFEA_error_text();
                    printf("\t\tError: In spring block, '%s' is not a recognised lvalue\n", lvalue);
                    printf("\t\tRecognised lvalues are:\n");
                    printf("\t\tspring_fname\n\n");
                    return FFEA_ERROR;
                }
            }
        } else {
            FFEA_error_text();
            printf("\t\tError: In system, '%s' is not a recognised tag\n", buf);
            printf("\t\tRecognised tags are:\n");
            printf("\t\tblob\n\tspring\n");
            return FFEA_ERROR;
        }
    }

    // Check we have correct number of blobs
    if (i != num_blobs) {
        FFEA_error_text();
        printf("Number of blobs read, %d, in is not the same as number specified in the script, %d\n", i, num_blobs);
        return FFEA_ERROR;
    }

    // And springs
    if (k > 1) {
        FFEA_error_text();
        printf("Read %d spring files. Only need 1, containing all possible springs!\n", k);
        return FFEA_ERROR;
    }

    // Add springs to world
    if (k == 1) {
        if (load_springs(spring_filename) != 0) {
            FFEA_error_text();
            printf("Problem building springs from %s\n", spring_filename);
            return FFEA_ERROR;
        }
    }

    return FFEA_OK;
}

/* */
void World::get_system_CoM(vector3 *system_CoM) {
    system_CoM->x = 0;
    system_CoM->y = 0;
    system_CoM->z = 0;
    scalar total_mass = 0;
    for (int i = 0; i < params.num_blobs; i++) {
        vector3 com;
        active_blob_array[i]->get_CoM(&com);
        system_CoM->x += com.x * active_blob_array[i]->get_mass();
        system_CoM->y += com.y * active_blob_array[i]->get_mass();
        system_CoM->z += com.z * active_blob_array[i]->get_mass();

        total_mass += active_blob_array[i]->get_mass();
    }
    system_CoM->x /= total_mass;
    system_CoM->y /= total_mass;
    system_CoM->z /= total_mass;
}

int World::load_springs(char *fname) {

    int i;
    FILE *in = NULL;
    const int max_line_size = 50;
    char line[max_line_size];

    // open the spring file
    if ((in = fopen(fname, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(fname)
    }
    printf("\t\tReading in springs file: %s\n", fname);

    // first line should be the file type "ffea springs file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of spring file\n")
    }
    if (strcmp(line, "ffea spring file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea spring file' (read '%s') \n", line)
    }

    // read in the number of springs in the file
    if (fscanf(in, "num_springs %d\n", &num_springs) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of springs\n")
    }
    printf("\t\t\tNumber of springs = %d\n", num_springs);

    // Allocate memory for springs
    spring_array = new Spring[num_springs];

    for (i = 0; i < num_springs; ++i) {
        if (fscanf(in, "%d %d %d %d %d %d %lf %lf\n", &spring_array[i].blob_index[0], &spring_array[i].conformation_index[0], &spring_array[i].node_index[0], &spring_array[i].blob_index[1], &spring_array[i].conformation_index[1], &spring_array[i].node_index[1], &spring_array[i].k, &spring_array[i].l) != 8) {
            FFEA_error_text();
            printf("Problem reading spring data from %s. Format is:\n\n", fname);
            printf("ffea spring file\nnum_springs ?\n");
            printf("blob_index_0 conformation_index_0 node_index_0 blob_index_1 conformation_index_1 node_index_1 k l\n");
            return FFEA_ERROR;
        }

    }

    fclose(in);
    printf("\t\t\tRead %d springs from %s\n", i, fname);
    activate_springs();
    return 0;
}

void World::activate_springs() {
    for (int i = 0; i < num_springs; ++i) {
        if (spring_array[i].conformation_index[0] == active_blob_array[spring_array[i].blob_index[0]]->conformation_index && spring_array[i].conformation_index[1] == active_blob_array[spring_array[i].blob_index[1]]->conformation_index) {
            spring_array[i].am_i_active = true;
        }
    }
}

void World::apply_springs() {
    scalar force_mag;
    vector3 n1, n0, force0, force1, sep, sep_norm;
    for (int i = 0; i < num_springs; ++i) {
        if (spring_array[i].am_i_active == true) {
            n1 = active_blob_array[spring_array[i].blob_index[1]]->get_node(spring_array[i].node_index[1]);
            n0 = active_blob_array[spring_array[i].blob_index[0]]->get_node(spring_array[i].node_index[0]);
            sep.x = n1.x - n0.x;
            sep.y = n1.y - n0.y;
            sep.z = n1.z - n0.z;
            sep_norm = normalise(&sep);
            force_mag = spring_array[i].k * (mag(&sep) - spring_array[i].l);

            force0.x = force_mag * sep_norm.x;
            force0.y = force_mag * sep_norm.y;
            force0.z = force_mag * sep_norm.z;

            force1.x = -1 * force_mag * sep_norm.x;
            force1.y = -1 * force_mag * sep_norm.y;
            force1.z = -1 * force_mag * sep_norm.z;

            active_blob_array[spring_array[i].blob_index[0]]->add_force_to_node(force0, spring_array[i].node_index[0]);
            active_blob_array[spring_array[i].blob_index[1]]->add_force_to_node(force1, spring_array[i].node_index[1]);
        }
    }
    return;
}

int World::get_next_script_tag(FILE *in, char *buf) {
    if (fscanf(in, "%*[^<]") != 0) {
        printf("White space removal error in get_next_script_tag(). Something odd has happened: this error should never occur...\n");
    }
    if (fscanf(in, "<%255[^>]>", buf) != 1) {
        FFEA_error_text();
        printf("Error reading tag in script file.\n");
        return FFEA_ERROR;
    }

    return FFEA_OK;
}

void World::apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N) {
    int i, j;
    for (i = 0; i < N; i++) {
        y[i] = 0;
        for (j = 0; j < N; j++)
            y[i] += M[i * N + j] * x[j];
    }
}

void World::do_es() {
    printf("Building BEM matrices\n");
    PB_solver.build_BEM_matrices();


    //	PB_solver.print_matrices();

    // Build the poisson matrices for each blob
    printf("Building Poisson matrices\n");
    for (int i = 0; i < params.num_blobs; i++) {
        active_blob_array[i]->build_poisson_matrices();
    }

    printf("Solving\n");
    for (int dual = 0; dual < 30; dual++) {


        // Perform Poisson solve step (K phi + E J_Gamma = rho)
        // Obtain the resulting J_Gamma
        int master_index = 0;
        for (int i = 0; i < params.num_blobs; i++) {
            active_blob_array[i]->solve_poisson(&phi_Gamma[master_index], &J_Gamma[master_index]);
            master_index += active_blob_array[i]->get_num_faces();
        }


        // Apply -D matrix to J_Gamma vector (work_vec = -D J_gamma)
        PB_solver.get_D()->apply(J_Gamma, work_vec);
        for (int i = 0; i < total_num_surface_faces; i++) {
            work_vec[i] *= -1;
        }

        // Solve for C matrix (C phi_Gamma = work_vec)
        nonsymmetric_solver.solve(PB_solver.get_C(), phi_Gamma, work_vec);

        scalar sum = 0.0, sumj = 0.0;
        for (int i = 0; i < total_num_surface_faces; i++) {
            sum += phi_Gamma[i];
            sumj += J_Gamma[i];
        }
        printf("\n");

        printf("<yophi> = %e <J> = %e\n", sum / total_num_surface_faces, sumj / total_num_surface_faces);

        //		for(i = 0; i < total_num_surface_faces; i++) {
        //			printf("<JPBSN> %d %e\n", dual, J_Gamma[i]);
        //		}

    }

    //	scalar sum = 0.0, sumj = 0.0;
    //	for(i = 0; i < total_num_surface_faces; i++) {
    //		sum += phi_Gamma[i];
    //		sumj += J_Gamma[i];
    //	}
    //	printf("<phi> = %e <J> = %e\n", sum/total_num_surface_faces, sumj/total_num_surface_faces);
    //
    //	for(i = 0; i < total_num_surface_faces; i++) {
    //		printf("<J> %e\n", J_Gamma[i]);
    //	}

    //	blob_array[0].print_phi();
}

void World::print_trajectory_and_measurement_files(int step, double wtime) {
    if ((step - 1) % (params.check * 10) != 0) {
        printf("step = %d\n", step);
    } else {
        printf("step = %d (simulation time = %.2es, wall clock time = %.2e hrs)\n", step, (scalar) step * params.dt, (omp_get_wtime() - wtime) / 3600.0);
    }

    vector3 system_CoM;
    get_system_CoM(&system_CoM);

    // Write traj and meas data
    if (measurement_out[params.num_blobs] != NULL) {
        fprintf(measurement_out[params.num_blobs], "%d\t", step);
    }
    if (stress_out != NULL) {
	fprintf(stderr, "Write to stress");
        fprintf(stress_out, "step\t%d\n", step);
    }

    for (int i = 0; i < params.num_blobs; i++) {

        // Write the node data for this blob
        fprintf(trajectory_out, "Blob %d, Conformation %d, step %d\n", i, active_conformation_index[i], step);
        active_blob_array[i]->write_nodes_to_file(trajectory_out);

        // Write the measurement data for this blob
        active_blob_array[i]->make_measurements(measurement_out[i], step, &system_CoM);

        // Output stress information
        active_blob_array[i]->make_stress_measurements(stress_out, i);

        // Output interblob_vdw info
        for (int j = i + 1; j < params.num_blobs; ++j) {
            active_blob_array[i]->calculate_vdw_bb_interaction_with_another_blob(measurement_out[params.num_blobs], j);
        }
    }

    fprintf(measurement_out[num_blobs], "\n");

    // Mark completed end of step with an asterisk (so that the restart code will know if this is a fully written step or if it was cut off half way through due to interrupt)
    fprintf(trajectory_out, "*\n");

    // Force all output in buffers to be written to the output files now
    fflush(trajectory_out);
    for (int i = 0; i < num_blobs + 1; ++i) {
        fflush(measurement_out[i]);
    }
}

void World::print_static_trajectory(int step, double wtime, int blob_index) {
    printf("Printing single trajectory of Blob %d for viewer\n", blob_index);
    // Write the node data for this blob
    fprintf(trajectory_out, "Blob %d, step %d\n", blob_index, step);
    active_blob_array[blob_index]->write_nodes_to_file(trajectory_out);
}
