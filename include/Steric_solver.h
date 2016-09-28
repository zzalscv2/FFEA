#include "VdW_solver.h"
#include <set>

class Steric_solver: public VdW_solver {
#ifdef USE_MPI
public:
  void do_interaction_mpi(Face *f1, Face *f2); 
#endif

private:
  void do_interaction(Face *f1, Face *f2); 

};
