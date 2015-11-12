#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cctype>
#include <omp.h>

#include <boost/program_options.hpp>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
//#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"

#define MAX_FNAME_LENGTH 255

using namespace std;
namespace b_po = boost::program_options;

int main(int argc, char *argv[])
{
	cout << "\n\n\n***************************************************\n\tFLUCTUATING FINITE ELEMENT ANALYSIS\n***************************************************\n\n" << endl;
	cout << " Version:\t" << FFEA_VERSION << " [" << FFEA_MASCOT << "]" << endl;
	cout << "Compiled:\t" << __DATE__ " at " << __TIME__ << endl;
	cout << "  Coding:\tRobin Richardson (pyrar@leeds.ac.uk), Ben Hanson (py09bh@leeds.ac.uk)\n" << endl;
	cout << "  Theory:\tOliver Harlen, Sarah Harris, Robin Oliver, Daniel Read, Robin Richardson, Ben Hanson\n" << endl;

	#ifdef FFEA_PARALLEL_WITHIN_BLOB
		cout << "Parallelisation switch: FFEA_PARALLEL_WITHIN_BLOB\n" << endl;
	#endif

	#ifdef FFEA_PARALLEL_PER_BLOB
		 cout << "Parallelisation switch: FFEA_PARALLEL_PER_BLOB\n" << endl;
	#endif
	
	// Return variable
	int myreturn = 0;
	
	// Get required args
	string script_fname;
	int mode;
	int frames_to_delete = 0;
	b_po::options_description desc("Allowed options");
	b_po::positional_options_description pdesc;
	b_po::variables_map var_map;
	try {

		// Options for visible and non-visible cmd line params
		desc.add_options()
			("help,h", "print usage message")
			("input-file,i", b_po::value<string>(&script_fname), "input script fname")
			("mode,m", b_po::value<int>(&mode)->default_value(0), "simulation mode (0 - FFEA, 1 - ENM, 2 - DMM, 3 - Timestep Calculator)")
			("delete_frames,l", b_po::value<int>(&frames_to_delete), "If restarting a simulation, this will delete the final 'arg' frames before restarting")
		;
		
		// 1 input file max! Option invisible (positional)
		pdesc.add("input-file", 1);

		// Parse
		b_po::store(b_po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), var_map);
		b_po::notify(var_map);

	// User related problems end here
	} catch(exception& e) {
		cerr << e.what() << endl;
	}

	// Boost makes sure options are valid. We check arguments are valid
	// --mode
	if(var_map.count("mode")) {
		if(mode != 0 && mode != 1 && mode != 2 && mode != 3 && mode != 4) {
			FFEA_error_text();		
			cout << "Unrecognised argument to --mode" << endl;
			cout << desc << endl;
			return -1;
		}
	}

	// --delete_frames
	if(var_map.count("delete_frames")) {
		if(frames_to_delete < 0) {
			FFEA_error_text();		
			cout << "Argument to --delete_frames cannot be less than zero. We can't magic up simulation frames for you! :D" << endl;
			cout << desc << endl;
			return -1;
		} else if (frames_to_delete == 0) {
			FFEA_error_text();		
			cout << "Argument to --delete_frames shouldn't equal zero as this won't have any effect and your simulation will crash again." << endl;
			cout << desc << endl;
			return -1;
		}
	}

	/* Starting actual functionality */
	// Help text is built in to boost	
	if (var_map.count("help")) {  
		cout << desc << endl;
		return 0;
	}

	// Check we have an input script
	if (var_map.count("input-file")) {  
		cout << "Input FFEA script - " << script_fname << "\n";
	} else {
		cout << "\n\nUsage: ffea [FFEA SCRIPT FILE (.ffea)] [OPTIONS]\n\n\n" << endl;
		cout << desc << endl;
		return FFEA_ERROR;
	}

	// The system of all proteins, electrostatics and water
	World *world;

	// Allocate the world
	world = new World();

	// Initialise the world, loading all blobs, parameters, electrostatics, etc.
	cout << "Initialising the world:\n" << endl;
	if(world->init(script_fname, frames_to_delete, mode) == FFEA_ERROR) {
		FFEA_error_text();
		cout << "Errors during initialisation mean World cannot be constructed properly." << endl;
		myreturn = FFEA_ERROR;

		// Delete the world (oh no!)
		cout << "Deleting world..." << endl;
		delete world;
		cout << "...done. World has been sucessfully destroyed." << endl;

		return myreturn;
	}
	
	// World is initialised. How shall we run FFEA?

	/* Calculate the maximum allowed timestep for this system */
	if(var_map.count("timestep")) {

		cout << endl << "FFEA - Calculating System Timescales" << endl << endl;
		if(world->get_smallest_time_constants() == FFEA_ERROR) {
			FFEA_error_text();
			cout << "Error occurred in 'get_smallest_time_constants()' in World\n" << endl;
			myreturn = FFEA_ERROR;
		} else {
			cout << "...done\n" << endl;
			myreturn = FFEA_OK;
		}

	}

	else if(mode == 1 || mode == 2 || mode == 3 || mode == 4) {
		
		if(mode == 1) {

			/* Elastic Network Model */
			cout << "\n\n\n***************************************************\n\tFFEA - Elastic Network Model\n***************************************************\n\n";
		} else if (mode == 2) {

			/* Dynamic Mode Model */
			cout << "\n\n\n***************************************************\n\tFFEA - Dynamic Mode Model\n***************************************************\n\n";
		} else if (mode == 4) {

			/* Rotne Prager Dynamic Mode Model */
			cout << "\n\n\n***************************************************\n\tFFEA - Rotne Prager Dynamic Mode Model\n***************************************************\n\n";
		}

		// Get the desired blobs
		cout << "\tFirstly, lets decide which blobs to analyse!" << endl << endl;
		set<int> blobs;
		int ablob;
		int error;
		char buf[5];
		while(true) {
			error = 0;
			try{
				cout << "\t\tEnter an index for the blob you would like an elastic network model for, or type 'q' to finish?:";
				scanf("%S", &buf);
				if(strcmp(buf, "q") == 0 or strcmp(buf, "Q") == 0) {
					cout << endl << "\tThat's all the blobs!" << endl;
					break;
				}

				for(int i = 0; i < 1; ++i) {
					if(isalpha(buf[i])) {
						FFEA_error_text();
						cout << "\tPlease enter a valid blob index (0 <= x <" << world->get_num_blobs() << ")" << endl;
						error = 1;
						break;
					}
				}
				
				if(error == 1) {
					continue;
				}

				ablob = atoi(buf);
				if(ablob < 0 || ablob >= world->get_num_blobs()) {
					FFEA_error_text();
					cout << "\tPlease enter a valid blob index (0 <= x <" << world->get_num_blobs() << ")" << endl;
					continue;
				} else {
					blobs.insert(ablob);
				}
			} catch (exception &e) {
				FFEA_error_text();
				cout << "Exception occured:" << endl;
				cout << e.what() << endl;
			}
		}

		if(blobs.size() == 0) {
			cout << "No blobs selected. No elastic network models will be found :(" << endl;
			return FFEA_OK;
		}

		// Get number of modes
		int num_modes;
		while(true) {
			try{
				cout << "\n\tHow many modes would you like to visualise?:";
				scanf("%s", &buf);
				for(int i = 0; i < 1; ++i) {
					if(isalpha(buf[i])) {
						FFEA_error_text();
						cout << "\tYou must choose a whole number i > 0 and i < 3N" << endl;
						error = 1;
						break;
					}
				}

				if(error == 1) {
					continue;
				}

				num_modes = atoi(buf);
				if(num_modes <= 0) {
					FFEA_error_text();
					cout << "\tYou must choose a whole number i > 0 and i < 3N" << endl;
					continue;
				}
				break;
			} catch(exception &e) {
				FFEA_error_text();
				cout << "Exception occured:" << endl;
				cout << e.what() << endl;
				cout << "\tPlease enter a whole number less than 3 * the total number of nodes in your system and greater than zero." << endl; 
			}
		}
		
		set<int>::iterator it;
		cout << endl << "\tYou've asked for " << num_modes << " modes to be calculated for blob(s):" << endl << "\t";
		for(it = blobs.begin(); it != blobs.end(); ++it) {
			cout << *it << " ";
		}

		if(mode == 1) {
			cout << endl << endl << "\tBeginning the calculation of the ENMs..." << endl;
			if(world->enm(blobs, num_modes) == FFEA_ERROR) {
				cout << endl;
				FFEA_error_text();
				cout << "Problem when calculating elastic network models." << endl;
			}
		} else if (mode == 2) {
			cout << endl << endl << "\tBeginning the calculation of the DMMs..." << endl;
			if(world->dmm(blobs, num_modes) == FFEA_ERROR) {
				cout << endl;
				FFEA_error_text();
				cout << "Problem when calculating dynamic mode models." << endl;
			}
		} else if (mode == 4) {
			cout << endl << endl << "\tBeginning the calculation of the Rotne-Prager DMMs..." << endl;
			if(world->dmm_rp(blobs, num_modes) == FFEA_ERROR) {
				cout << endl;
				FFEA_error_text();
				cout << "Problem when calculating RP dynamic mode models." << endl;
			}
		}
		cout << "\tdone!" << endl;

	} else if(mode == 0) {
		
		/* Full FFEA */

		// Run the world for the specified number of time steps
		cout << "Running world...\n" << endl;
		if(world->run() == FFEA_ERROR) {
			FFEA_error_text();
			cout << "Error occurred when running World\n" << endl;
			myreturn = FFEA_ERROR;
		} else {
			cout << "...done\n" << endl;
			myreturn = FFEA_OK;
		}
	} else {
		cout << "Error. Unrecognised argument to --mode \"" << mode << "\"" << endl << endl;
		cout << desc << endl;
	}

	// Delete the world (oh no!)
	cout << "Deleting world..." << endl;
	delete world;
	cout << "...done. World has been sucessfully destroyed." << endl;

	return myreturn;
}
