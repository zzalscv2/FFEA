#include "VdW_solver.h"
#include <set>

class Steric_solverII: public VdW_solver {

private:
  void do_interaction(Face *f1, Face *f2); 

};
