#include "VdW_solver.h"
#include <set>

class Steric_solver: public VdW_solver {

private:
  void do_interaction(Face *f1, Face *f2);

  void do_interaction(Face *f1, Face *f2, scalar *blob_corr);

};
