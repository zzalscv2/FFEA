#include <cstdio>
#include <cmath>

using namespace std;

int main(int argc char **argv) {

	// Check input parameters
	if(argc != 4) {
		cout << "Usage: ./surface_coarse_grainer_12122014 [INPUT .surf FNAME] [OUTPUT FNAME] [COARSENESS LEVEL (angstroms)]" << endl;
		return -1;
	}

	// Create and initialise surface
	Surface surf;
	surf.init(argv[1]);
	surf.coarsen(argv[3]);
	surf.write_to_file(argv[2]);
	return 0;
}
