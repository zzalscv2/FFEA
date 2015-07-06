#include <iostream>
#include <stdlib.h>
#include <ctime>

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
	b_po::options_description desc("Allowed options");
	b_po::positional_options_description pdesc;
	b_po::variables_map var_map;
	try {

		// Options for visible and non-visible cmd line params
		desc.add_options()
			("help,h", "print usage message")
			("input-file,i", b_po::value(&script_fname), "input script fname")
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

	// Help text is built in to boost	
	if (var_map.count("help")) {  
		cout << desc << "\n";
		return 0;
	}

	// Check we have an input script
	if (var_map.count("input-file")) {  
		cout << "Input FFEA script - " << var_map["input-file"].as<string>() << "\n";
	} else {
		cout << "\n\nUsage: ffea [FFEA SCRIPT FILE (.ffea)] [OPTIONS]\n\n\n" << endl;
		return FFEA_ERROR;
	}

	// The system of all proteins, electrostatics and water
	World *world;

	// Allocate the world
	world = new World();

	// Initialise the world, loading all blobs, parameters, electrostatics, etc.
	cout << "Initialising the world:\n" << endl;
	if(world->init(script_fname) == FFEA_ERROR) {
		FFEA_error_text();
		cout << "Errors during initialisation mean World cannot be constructed properly." << endl;
		myreturn = FFEA_ERROR;
	} else {

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
	}

	// Delete the world (oh no!)
	cout << "Deleting world..." << endl;
	delete world;
	cout << "...done. World has been sucessfully destroyed." << endl;

	return myreturn;
}
