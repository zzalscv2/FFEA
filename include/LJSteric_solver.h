#include "VdW_solver.h"
#include <set>

class LJSteric_solver: public VdW_solver {

private:
  void do_interaction(Face *f1, Face *f2); 

};
