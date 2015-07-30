#include "PreComp_solver.h"

using namespace std;

PreComp_solver::PreComp_solver() { 
  msgc = 0;
} 

/** @brief destrctor: deallocates pointers. */
PreComp_solver::~PreComp_solver() {
  delete U;
  delete F;
}


int PreComp_solver::msg(int whatever){
  cout << "--- PreComp " << msgc << ": " << whatever << endl;
  msgc += 1;
  return 0; 
} 

int PreComp_solver::msg(string whatever){
  cout << "--- PreComp " << msgc << ": " << whatever << endl;
  msgc += 1;
  return 0; 
} 


/** 
 * @brief Read input precomputed tables
 * @param[in] vector<string> types: types of beads present.
 * @param[in] int inputData: 1 means read .force and .pot files,
 *                 while 2 means read .pot and calculate the forces
 * @param[in] int approach: if .vdw -> use .vdw files and solve faces;
 *                          if .solid -> use .solid files and solve volumes.
 * @detail read <type_i>-<type_j>.force and <type_i>-<type_j>.pot files
 *         for all the possible pairs, setting up the potential functions.
 *         All the .pot and .force files need to have the same x_range
 *         same Dx, and same number of points. It is allowed, however, that 
 *         they have a different number of lines at the beginning starting with "#". 
 */        
int PreComp_solver::init(PreComp_params &pc_params) {

   /* Firstly, we get the number of lines, x_range, 
    * and Dx from the first pair type potential file. 
    * We check that the x_range and Dx 
    * are the same for all the files. Secondly, 
    * we allocate the F and U arrays, and then
    * we store the y_values, re-reading the files.
    * Finally, the functions get_U and get_F will be ready. 
    */  
   
   if (pc_params.types.size() == 0) return 0;

   stringstream ssfile; 
   ifstream fin;
   string line; 
   vector<string> vec_line; 
   scalar x_0, x;
   
   /*--------- FIRSTLY ----------*/ 
   /* get the number of interactions */
   for (int i=1; i<pc_params.types.size() + 1; i++){
     nint += i;
   } 
   msg(nint);
   // and the number of bead types: 
   ntypes = pc_params.types.size();

   /* read the first file, and get the number of lines
    * (n_values), Dx, and x_range */
   // the first file:
   ssfile << pc_params.folder << "/" << pc_params.types[0] << "-" << pc_params.types[0] << ".pot";
   fin.open(ssfile.str());
   // get the first line that does not start with "#"
   getline(fin, line);
   while (line.find("#", 0, 1) == 0) {
     getline(fin, line); 
   } 
   // get x_0, i. e., parse the line:
   boost::split( vec_line, line, boost::is_any_of(" \t"));
   x_0 = stod(vec_line[0]); 
   // and store it: 
   x_range[0] = x_0;
   getline(fin, line);
   boost::split( vec_line, line, boost::is_any_of(" \t"));
   x = stod(vec_line[0]); 
   Dx = x - x_0;
   n_values = 2;
   x_0 = x;
   while ( getline(fin, line) ) {
      n_values += 1;
      boost::split( vec_line, line, boost::is_any_of(" \t"));
      x = stod(vec_line[0]); 
      if (fabs(Dx - x + x_0) > 1e-6) { 
        FFEA_error_text();
        msg("Aborting: delta x was found to be non-uniform in:");
        msg(ssfile.str());
        return FFEA_ERROR;
      } 
      x_0 = x; 
   }  
   fin.close();

   /*--------- SECONDLY ----------*/ 
   // allocate:
   U = new scalar[n_values * nint];    
   F = new scalar[n_values * nint];    
      
   // and load potentials and forces:
   read_tabulated_values(pc_params, "pot", U);
   if (pc_params.inputData == 1) {
     read_tabulated_values(pc_params, "force", F);
   } else if (pc_params.inputData == 2) {
     calc_force_from_pot();
   } else {
     FFEA_error_text();
     msg("invalid value for precomp->inputData");
     return FFEA_ERROR;
   }


   return FFEA_OK; 
}


/** Calculate F as the numerical derivative of U, instead of reading .force files 
  * @param[in] int total_int: total number of different interactions.
  */
int PreComp_solver::calc_force_from_pot() {

   scalar ym, y0, yM;
   scalar twoDx = 2*Dx;

   int index = -1;
   for (int i=0; i<nint; i++) {
     index += 1; 
     F[index] = (U[index] - U[index + 1]) / Dx;
     for (int j=1; j<n_values-1; j++) {
       index += 1;
       F[index] = (U[index -1] - U[index + 1]) / twoDx; 
     }
     index += 1;
     F[index] = 0; 
   } 
   return FFEA_OK;


}

/** Read either the .pot or .force files
  * and store its contents either in U or F */
int PreComp_solver::read_tabulated_values(PreComp_params &pc_params, string kind, scalar *Z){

   stringstream ssfile; 
   ifstream fin;
   string line; 
   vector<string> vec_line; 
   scalar x_0, x;

   int m_values;
   int index = 0;
   // read the potentials/forces in while checking:
   for (int i=0; i<pc_params.types.size(); i++) {
     for (int j=i; j<pc_params.types.size(); j++) { 
       ssfile << pc_params.folder << "/" << pc_params.types[i] << "-" << pc_params.types[j] << "." << kind;
       fin.open(ssfile.str());
       // get the first line that does not start with "#"
       getline(fin, line);
       while (line.find("#", 0, 1) == 0) {
         getline(fin, line); 
       }  
       // get x_0, i. e., parse the line:
       boost::split( vec_line, line, boost::is_any_of(" \t"));
       x_0 = stod(vec_line[0]); 
       // check it: 
       if (x_range[0] != x_0){
          FFEA_error_text();
          msg("different starting point for file:");
          msg(ssfile.str());
          return FFEA_ERROR;
       }
       // and store the potential:
       Z[index] = stod(vec_line[1]);
       index += 1;
       // get the next line, and check Dx:
       getline(fin, line);
       boost::split( vec_line, line, boost::is_any_of(" \t"));
       x = stod(vec_line[0]); 
       if (fabs(Dx - x + x_0) > 1e-6){
          FFEA_error_text();
          msg("different step size for x in file:");
          msg(ssfile.str());
          return FFEA_ERROR;
       }
       // and store the potential:
       Z[index] = stod(vec_line[1]);
       index += 1;
       // Now, go for the rest of the file:
       m_values = 2;
       x_0 = x;
       while ( getline(fin, line) ) {
          m_values += 1;
          boost::split( vec_line, line, boost::is_any_of(" \t"));
          x = stod(vec_line[0]); 
          // check Dx at every line:
          if (fabs(Dx - x + x_0) > 1e-6) { 
            FFEA_error_text();
            msg("Aborting; delta x was found to be non-uniform for file:");
            msg(ssfile.str());
            return FFEA_ERROR;
          } 
          x_0 = x; 
          // and check that it is not too long:
          if (m_values > n_values) { 
            FFEA_error_text();
            msg("Aborting; too many points for file:");
            msg(ssfile.str());
            return FFEA_ERROR;
          } 
          Z[index] = stod(vec_line[1]);
          index += 1;
       }  
       // and check that it is not too short:
       if (m_values != n_values) { 
         FFEA_error_text();
         msg("Aborting; too many points for file:");
         msg(ssfile.str());
         return FFEA_ERROR;
       } 
       fin.close();
       ssfile.str("");
     }
   } 
     
   return 0;
}

/* @detail Compute (through finterpolate) the value of the potential at x 
 * between bead types typei and typej. Kept just for clarity, may be removed
 * in the future.
 */
scalar PreComp_solver::get_U(scalar x, int typei, int typej) {
 
  return finterpolate(U, x, typei, typej); 

}

/* @detail Compute (through finterpolate) the value of the potential at x 
 * between bead types typei and typej. Kept just for clarity, may be removed
 * in the future.
 */
scalar PreComp_solver::get_F(scalar x, int typei, int typej) {
 
  return finterpolate(F, x, typei, typej); 

}

/** @brief Get the value Z value (either U or F) at x between types typei and typej.
  * @detail The function is currently private and accessed through get_U and get_F.
  * For production purposes it may be reasonable to move it to public, and 
  * remove get_U and get_F
  */
scalar PreComp_solver::finterpolate(scalar *Z, scalar x, int typei, int typej){

   int index;
   scalar y0, y1, x0, x1;

   // sort so that typei <= typej
   if (typei > typej){
     index = typei;
     typei = typej;
     typej = index;
   } 


   // get the index for the closest (bottom) value:
   index = typej - typei;
   for (int i=ntypes; i>typei; i--) {
      index += i;
   }
   index = index*n_values;
   index += x/Dx ;

   // interpolate:
   x0 = int(x/Dx) * Dx; 
   x1 = x0 + Dx; 
   msg("");  
   return Z[index] + (Z[index +1] - Z[index])*(x - x0)/Dx;

}  
