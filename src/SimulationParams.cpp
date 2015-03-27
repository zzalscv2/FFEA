#include "SimulationParams.h"


		SimulationParams::SimulationParams()
		{
			dt = 0;
			num_steps = -1;
			check = 0;
			num_blobs = 0;
			num_conformations = NULL;
			num_states = NULL;
			rng_seed = 0;
			kT = 0;
			max_iterations_cg = 0;
			epsilon2 = 0;
			es_update = 0;
			es_N_x = -1;
			es_N_y = -1;
			es_N_z = -1;
			kappa = 0;
			dielec_ext = 0;
			epsilon_0 = 0;
			restart = 0;
			calc_vdw = -1;
			calc_es = 0;
			calc_noise = 1;

			wall_x_1 = WALL_TYPE_PBC;
			wall_x_2 = WALL_TYPE_PBC;
			wall_y_1 = WALL_TYPE_PBC;
			wall_y_2 = WALL_TYPE_PBC;
			wall_z_1 = WALL_TYPE_PBC;
			wall_z_2 = WALL_TYPE_PBC;

			sticky_wall_xz = 0;

			do_stokes = 0;
			stokes_visc = -1;

			vdw_r_eq = -1;
			vdw_eps = -1;

			trajectory_out_fname_set = 0;
			measurement_out_fname_set = 0;
			stress_out_fname_set = 0;

			vdw_params_fname_set = 0;

			sprintf(trajectory_out_fname, "\n");
			sprintf(temp_fname, "\n");
			sprintf(stress_out_fname, "\n");
			measurement_out_fname = NULL;
			sprintf(vdw_params_fname, "\n");
		}


		SimulationParams::~SimulationParams()
		{
			dt = 0;
			num_steps = -1;
			check = 0;
			num_blobs = 0;
			delete[] num_conformations;
			num_conformations = NULL;
			delete[] num_states;
			num_states = NULL;
			rng_seed = 0;
			kT = 0;
			max_iterations_cg = 0;
			epsilon2 = 0;
			es_update = 0;
			es_N_x = -1;
			es_N_y = -1;
			es_N_z = -1;
			kappa = 0;
			dielec_ext = 0;
			epsilon_0 = 0;
			restart = 0;
			calc_vdw = -1;
			calc_es = 0;
			calc_noise = 0;

			wall_x_1 = -1;
			wall_x_2 = -1;
			wall_y_1 = -1;
			wall_y_2 = -1;
			wall_z_1 = -1;
			wall_z_2 = -1;

			sticky_wall_xz = 0;

			do_stokes = 0;
			stokes_visc = -1;

			vdw_r_eq = -1;
			vdw_eps = -1;

			trajectory_out_fname_set = 0;
			measurement_out_fname_set = 0;
			stress_out_fname_set = 0;
			vdw_params_fname_set = 0;

			sprintf(trajectory_out_fname, "\n");
			sprintf(temp_fname, "\n");
			delete[] measurement_out_fname;
			sprintf(stress_out_fname, "\n");
			sprintf(vdw_params_fname, "\n");
		}

		int SimulationParams::validate()
		{
			printf("Validating parameters...\n");

			if(restart != 0 && restart != 1) {
				FFEA_ERROR_MESSG("Required: Restart flag, 'restart', must be 0 (false) or 1 (true).\n");
			}

			if(num_steps < 0) {
				FFEA_ERROR_MESSG("Required: Number of time steps, 'num_steps', must be greater than or equal to 0.\n");
			}
			if(num_blobs <= 0) {
				FFEA_ERROR_MESSG("\tRequired: Number of Blobs, 'num_blobs', must be greater than 0.\n");
			}
			for(int i = 0; i < num_blobs; ++i) {
				if(num_conformations[i] <= 0) {
					FFEA_ERROR_MESSG("\tRequired: Number of Conformations, 'num_conformations[%d]', must be greater than 0.\n", i);
				}
				if(num_states[i] <= 0) {
					FFEA_ERROR_MESSG("\tRequired: Number of States, 'num_states[%d]', must be greater than 0.\n", i);
				}
				if(num_conformations[i] > num_states[i]) {
					FFEA_ERROR_MESSG("\tRequired: Number of Conformations, 'num_conformations[%d]', must be less than or equal to Number of States, 'num_states[%d]'.\n", i, i);
				}
			}
			if(es_N_x < 1) {
				printf("FRIENDLY WARNING: Length of the nearest neighbour lookup grid, 'es_N_x', is less than 1. No electrostatics will be simulated.\n");
			} else if(es_N_y < 1) {
				printf("FRIENDLY WARNING: Length of the nearest neighbour lookup grid, 'es_N_y', is less than 1. No electrostatics will be simulated.\n");
			} else if(es_N_z < 1) {
				printf("FRIENDLY WARNING: Length of the nearest neighbour lookup grid, 'es_N_z', is less than 1. No electrostatics will be simulated.\n");
			} else {
				if(es_h <= 0) {
					FFEA_ERROR_MESSG("Required: Nearest neighbour lookup grid cell dimension, 'es_h', must be greater than 0.\n");
				}
			}

			if(kappa < 0) {
				FFEA_ERROR_MESSG("Required: Inverse Debye Screening length, 'kappa', cannot be negative.\n");
			}

			if(dielec_ext <= 0) {
				FFEA_ERROR_MESSG("Required: Exterior dielectric constant, 'dielec_ext', must be greater than 0.\n");
			}

			if(epsilon_0 <= 0) {
				FFEA_ERROR_MESSG("Required: Permittivity of free space, 'epsilon_0', must be greater than 0.\n");
			}

			if(calc_vdw != 0 && calc_vdw != 1) {
				FFEA_ERROR_MESSG("Required: 'calc_vdw', must be 0 (no) or 1 (yes).\n");
			}

			if(calc_vdw == 1) {
				if(vdw_params_fname_set == 0) {
					FFEA_ERROR_MESSG("VdW forcefield params file name required (vdw_forcefield_params).\n");
				}
			}

			if(calc_es != 0 && calc_es != 1) {
				FFEA_ERROR_MESSG("Required: 'calc_es', must be 0 (no) or 1 (yes).\n");
			}

			if(calc_noise != 0 && calc_noise != 1) {
				FFEA_ERROR_MESSG("Required: 'calc_noise', must be 0 (no) or 1 (yes).\n");
			}

			if(trajectory_out_fname_set == 0) {
				FFEA_ERROR_MESSG("Trajectory output file name required (trajectory_out_fname).\n");
			}

			if(measurement_out_fname_set == 0) {
				FFEA_ERROR_MESSG("Measurement output file name required (measurement_out_fname).\n");
			} else {
				measurement_out_fname = new char*[num_blobs + 1];
				char temp[MAX_FNAME_SIZE];
				char base[MAX_FNAME_SIZE];
				char ext[MAX_FNAME_SIZE];
				char *pch = NULL;
				sprintf(temp, temp_fname);
				pch = strtok(temp, ".");
				int num_points = -1;
				while(pch != NULL) {
					num_points++;
					pch = strtok(NULL, ".");
				}
				sprintf(temp, temp_fname);
				pch = strtok(temp, ".");
				sprintf(base, pch);
				for(int i = 1; i < num_points; ++i) {
					pch = strtok(NULL, ".");
					strcat(base, pch);
				}

				sprintf(ext, ".");
				if(num_points == 0) {
					strcat(ext, "out");
				} else {
					pch = strtok(NULL, ".");
					strcat(ext, pch);
				}
				
				for(int i = 0; i < num_blobs + 1; ++i) {
					measurement_out_fname[i] = new char[MAX_FNAME_SIZE];
					sprintf(measurement_out_fname[i], base);
					if(i < num_blobs) {
						sprintf(temp, "_blob%d", i);
						strcat(measurement_out_fname[i], temp);
					} else {
						strcat(measurement_out_fname[i], "_world");	
					}
					strcat(measurement_out_fname[i], ext);
				}
			}

			if(do_stokes == 1 && stokes_visc <= 0) {
				FFEA_ERROR_MESSG("do_stokes flag is set, so stokes_visc must be set to a value greater than 0.\n");
			}

			printf("...done\n");

			printf("Parameters:\n");
			printf("\trestart = %d\n", restart);
			printf("\tdt = %e\n", dt);
			printf("\tnum_steps = %lld\n", num_steps);
			printf("\tcheck = %d\n", check);
			printf("\tnum_blobs = %d\n", num_blobs);
			for(int i = 0; i < num_blobs; ++i) {
				printf("\tBlob %d:\n", i);
				printf("\t\tnum_conformations = %d\n", num_conformations[i]);
				printf("\t\tnum_states = %d\n", num_states[i]);
			}
			printf("\trng_seed = %d\n", rng_seed);
			printf("\tkT = %e\n", kT);
			printf("\ttrajectory_out_fname = %s\n", trajectory_out_fname);
			for(int i = 0; i < num_blobs + 1; ++i) {
				printf("\tmeasurement_out_fname %d = %s\n", i, measurement_out_fname[i]);
			}
			printf("\tstress_out_fname = %s\n", stress_out_fname);
			printf("\tvdw_forcefield_params = %s\n", vdw_params_fname);
			printf("\tmax_iterations_cg = %d\n", max_iterations_cg);
			printf("\tepsilon2 = %e\n", epsilon2);
			printf("\tes_update = %d\n", es_update);
			printf("\tes_N_x = %d\n", es_N_x);
			printf("\tes_N_y = %d\n", es_N_y);
			printf("\tes_N_z = %d\n", es_N_z);
			printf("\tes_h = %e x inverse kappa\n", es_h);
			printf("\tkappa = %e\n", kappa);
			printf("\tepsilon_0 = %e\n", epsilon_0);
			printf("\tdielec_ext = %e\n", dielec_ext);
			printf("\tcalc_vdw = %d\n", calc_vdw);
			printf("\tcalc_es = %d\n", calc_es);
			printf("\tcalc_noise = %d\n", calc_noise);
			printf("\tdo_stokes = %d\n", do_stokes);
			printf("\tstokes_visc = %f\n", stokes_visc);

			return FFEA_OK;
		}

		/*
		 * Expects a string of the form "lvalue = rvalue" where lvalue must be a recognised parameter name.
		 */
		int SimulationParams::parse_param_assignment(char *str)
		{
			const int buf_size = 101;
			char lvalue[buf_size];
			char rvalue[buf_size];
			int rv;

			rv = sscanf(str, "%100[^=]=%s", lvalue, rvalue);

			if(rv != 2) {
				FFEA_ERROR_MESSG("\tError parsing parameter assignment, '%s'\n", str);
			}

			rv = sscanf(lvalue, "%s", lvalue);

			// Carry out parameter assignments
			if(strcmp(lvalue, "dt") == 0) {
				dt = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, dt);
			} else if(strcmp(lvalue, "epsilon") == 0) {
				epsilon2 = atof(rvalue);
				epsilon2 *= epsilon2;
				printf("\tSetting %s = %g\n", lvalue, atof(rvalue));
			} else if(strcmp(lvalue, "num_steps") == 0) {
				// convert to float first so that user may write numbers of the form 1e4 for 10000 etc
				num_steps = (long long)(atof(rvalue));
				printf("\tSetting %s = %lld\n", lvalue, num_steps);
			} else if(strcmp(lvalue, "max_iterations_cg") == 0) {
				max_iterations_cg = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, max_iterations_cg);
			} else if(strcmp(lvalue, "check") == 0) {
				check = (int)atof(rvalue);
				printf("\tSetting %s = %d\n", lvalue, check);
			} else if(strcmp(lvalue, "num_blobs") == 0) {
				num_blobs = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, num_blobs);
			} else if(strcmp(lvalue, "num_conformations") == 0) {
				if(num_blobs == 0) {
					FFEA_ERROR_MESSG("Error. num_blobs must be specified before num_conformations in the input scripts\n")
				}
				num_conformations = new int[num_blobs];
				char * pch;
				for(int i = 0; i < num_blobs; ++i) {
					if(i == 0) {
						pch = strtok(rvalue, "(,");
					} else {
						pch = strtok(NULL, " ,)");
					}	
					if(pch == NULL) {
						FFEA_ERROR_MESSG("num_conformations must have %d (num_blobs) arguments in the form (%%d,%%d,...,%%d). You've only specified %d arguments.\n", num_blobs, i);
					}
					num_conformations[i] = atoi(pch);
					printf("\tSetting Blob %d, %s = %d\n", i, lvalue, num_conformations[i]);
				}

				// Check that we are closing the bracket
				pch = strtok(NULL, " ,)");
				if(pch != NULL) {
					FFEA_ERROR_MESSG("num_conformations must have %d (num_blobs) arguments in the form (%%d,%%d,...,%%d). You've specified at least %d arguments.\n", num_blobs, num_blobs + 1);
				}

			} else if(strcmp(lvalue, "num_states") == 0) {
				if(num_blobs == 0) {
					FFEA_ERROR_MESSG("Error. num_blobs must be specified before num_states in the input scripts\n")
				}
				num_states = new int[num_blobs];
				char * pch;
				for(int i = 0; i < num_blobs; ++i) {
					if(i == 0) {
						pch = strtok(rvalue, "(,");
					} else {
						pch = strtok(NULL, " ,)");
					}
					if(pch == NULL) {
						FFEA_ERROR_MESSG("num_states must have %d (num_blobs) arguments in the form (%%d,%%d,...,%%d). You've only specified %d arguments.\n", num_blobs, i);
					}
					num_states[i] = atoi(pch);
					printf("\tSetting Blob %d, %s = %d\n", i, lvalue, num_states[i]);
				}
				
				// Check that we are closing the bracket
				pch = strtok(NULL, " ,)");
				if(pch != NULL) {
					FFEA_ERROR_MESSG("num_states must have %d (num_blobs) arguments in the form (%%d,%%d,...,%%d). You've specified at least %d arguments.\n", num_blobs, num_blobs + 1);
				}

			} else if(strcmp(lvalue, "es_update") == 0) {
				es_update = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, es_update);
			} else if(strcmp(lvalue, "es_N_x") == 0) {
				es_N_x = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, es_N_x);
			} else if(strcmp(lvalue, "es_N_y") == 0) {
				es_N_y = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, es_N_y);
			} else if(strcmp(lvalue, "es_N_z") == 0) {
				es_N_z = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, es_N_z);
			} else if(strcmp(lvalue, "sticky_wall_xz") == 0) {
				sticky_wall_xz = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, sticky_wall_xz);
			} else if(strcmp(lvalue, "es_h") == 0) {
				es_h = atof(rvalue);
				printf("\tSetting %s = %e x inverse kappa\n", lvalue, es_h);
			} else if(strcmp(lvalue, "rng_seed") == 0) {
				if(strcmp(rvalue, "time") == 0)
					rng_seed = time(NULL);
				else
					rng_seed = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, rng_seed);
			} else if(strcmp(lvalue, "kT") == 0) {
				kT = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, kT);
			} else if(strcmp(lvalue, "kappa") == 0) {
				kappa = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, kappa);
			} else if(strcmp(lvalue, "dielec_ext") == 0) {
				dielec_ext = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, dielec_ext);
			} else if(strcmp(lvalue, "epsilon_0") == 0) {
				epsilon_0 = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, epsilon_0);
			} else if(strcmp(lvalue, "restart") == 0) {
				restart = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, restart);
			} else if(strcmp(lvalue, "calc_vdw") == 0) {
				calc_vdw = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, calc_vdw);
			} else if(strcmp(lvalue, "calc_es") == 0) {
				calc_es = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, calc_es);
			} else if(strcmp(lvalue, "calc_noise") == 0) {
				calc_noise = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, calc_noise);
			} else if(strcmp(lvalue, "do_stokes") == 0) {
				do_stokes = atoi(rvalue);
				printf("\tSetting %s = %d\n", lvalue, do_stokes);
			} else if(strcmp(lvalue, "stokes_visc") == 0) {
				stokes_visc = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, stokes_visc);
			} else if(strcmp(lvalue, "wall_x_1") == 0) {
				if(strcmp(rvalue, "PBC") == 0) {
					wall_x_1 = WALL_TYPE_PBC;
				} else if(strcmp(rvalue, "HARD") == 0) {
					wall_x_1 = WALL_TYPE_HARD;
				} else if(strcmp(rvalue, "STOP") == 0) {
					wall_x_1 = WALL_TYPE_STOP;
				}
				printf("\tSetting %s = %d\n", lvalue, wall_x_1);
			} else if(strcmp(lvalue, "wall_x_2") == 0) {
				if(strcmp(rvalue, "PBC") == 0) {
					wall_x_2 = WALL_TYPE_PBC;
				} else if(strcmp(rvalue, "HARD") == 0) {
					wall_x_2 = WALL_TYPE_HARD;
				} else if(strcmp(rvalue, "STOP") == 0) {
					wall_x_2 = WALL_TYPE_STOP;
				}
				printf("\tSetting %s = %d\n", lvalue, wall_x_2);
			} else if(strcmp(lvalue, "wall_y_1") == 0) {
				if(strcmp(rvalue, "PBC") == 0) {
					wall_y_1 = WALL_TYPE_PBC;
				} else if(strcmp(rvalue, "HARD") == 0) {
					wall_y_1 = WALL_TYPE_HARD;
				} else if(strcmp(rvalue, "STOP") == 0) {
					wall_y_1 = WALL_TYPE_STOP;
				}
				printf("\tSetting %s = %d\n", lvalue, wall_y_1);
			} else if(strcmp(lvalue, "wall_y_2") == 0) {
				if(strcmp(rvalue, "PBC") == 0) {
					wall_y_2 = WALL_TYPE_PBC;
				} else if(strcmp(rvalue, "HARD") == 0) {
					wall_y_2 = WALL_TYPE_HARD;
				} else if(strcmp(rvalue, "STOP") == 0) {
					wall_y_2 = WALL_TYPE_STOP;
				}
				printf("\tSetting %s = %d\n", lvalue, wall_y_2);
			} else if(strcmp(lvalue, "wall_z_1") == 0) {
				if(strcmp(rvalue, "PBC") == 0) {
					wall_z_1 = WALL_TYPE_PBC;
				} else if(strcmp(rvalue, "HARD") == 0) {
					wall_z_1 = WALL_TYPE_HARD;
				} else if(strcmp(rvalue, "STOP") == 0) {
					wall_z_1 = WALL_TYPE_STOP;
				}
				printf("\tSetting %s = %d\n", lvalue, wall_z_1);
			} else if(strcmp(lvalue, "wall_z_2") == 0) {
				if(strcmp(rvalue, "PBC") == 0) {
					wall_z_2 = WALL_TYPE_PBC;
				} else if(strcmp(rvalue, "HARD") == 0) {
					wall_z_2 = WALL_TYPE_HARD;
				} else if(strcmp(rvalue, "STOP") == 0) {
					wall_z_2 = WALL_TYPE_STOP;
				}
				printf("\tSetting %s = %d\n", lvalue, wall_z_2);
			} else if(strcmp(lvalue, "vdw_r_eq") == 0) {
				vdw_r_eq = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, vdw_r_eq);
			} else if(strcmp(lvalue, "vdw_eps") == 0) {
				vdw_eps = atof(rvalue);
				printf("\tSetting %s = %e\n", lvalue, vdw_eps);
			} else if(strcmp(lvalue, "trajectory_out_fname") == 0) {
				if(strlen(rvalue) >= MAX_FNAME_SIZE) {
					FFEA_ERROR_MESSG("trajectory_out_fname is too long. Maximum filename length is %d characters.\n", MAX_FNAME_SIZE - 1)
				}
				sprintf(trajectory_out_fname, "%s", rvalue);
				trajectory_out_fname_set = 1;
				printf("\tSetting %s = %s\n", lvalue, trajectory_out_fname);
			} else if(strcmp(lvalue, "measurement_out_fname") == 0) {
				if(strlen(rvalue) >= MAX_FNAME_SIZE) {
					FFEA_ERROR_MESSG("measurement_out_fname is too long. Maximum filename length is %d characters.\n", MAX_FNAME_SIZE - 1)
				}
				sprintf(temp_fname, "%s", rvalue);
				measurement_out_fname_set = 1;
				printf("\tSetting %s = %s\n", lvalue, temp_fname);
			} else if(strcmp(lvalue, "stress_out_fname") == 0) {
				if(strlen(rvalue) >= MAX_FNAME_SIZE) {
					FFEA_ERROR_MESSG("stress_out_fname is too long. Maximum filename length is %d characters.\n", MAX_FNAME_SIZE - 1)
				}
				sprintf(stress_out_fname, "%s", rvalue);
				stress_out_fname_set = 1;
				printf("\tSetting %s = %s\n", lvalue, stress_out_fname);
			} else if(strcmp(lvalue, "vdw_forcefield_params") == 0) {
				if(strlen(vdw_params_fname) >= MAX_FNAME_SIZE) {
					FFEA_ERROR_MESSG("vdw_forcefield_params is too long. Maximum filename length is %d characters.\n", MAX_FNAME_SIZE - 1)
				}
				sprintf(vdw_params_fname, "%s", rvalue);
				vdw_params_fname_set = 1;
				printf("\tSetting %s = %s\n", lvalue, vdw_params_fname);	
			} else {
				FFEA_error_text();
				printf("\tError: '%s' is not a recognised lvalue\n", lvalue);
				printf("\tRecognised lvalues are:\n");
				printf("\tdt\n\tepsilon\n\tnum_steps\n\tmax_iterations_cg\n\tcheck\n\tes_update\n\ttrajectory_out_fname\n\tmeasurement_out_fname\n\tstress_out_fname\n\tes_N_x\n\tes_N_y\n\tes_N_z\n\tes_h\n\trng_seed\n\tkT\n\tkappa\n\tdielec_ext\tepsilon_0\n\trestart\n\tcalc_vdw\n\tcalc_noise\n\n");
				return FFEA_ERROR;
			}

			return FFEA_OK;
		}
