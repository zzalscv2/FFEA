// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

/*
 *      rod_structure.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk
 */

#include <iostream>
#include <fstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>

#include <stdio.h>
#include <random>
#include "rod_structure.h"

namespace rod {

/** Easy access to 1-d arrays */
#define odx(x) x*3
#define ody(y) (y*3)+1
#define odz(z) (z*3)+2
#define bend_index 0
#define stretch_index 1
#define twist_index 2

     /**---------**/
    /** Utility **/
   /**---------**/

/**
 Convert a vector containing strings to a vector of floats.
*/
std::vector<float> stof_vec(std::vector<std::string> vec_in){
    std::vector<float> vec_out(vec_in.size());
    for (unsigned int i=0; i<vec_in.size(); i++){
         vec_out[i] = std::stof(vec_in[i]);
    }
    return vec_out;
}

float random_number(float A, float B, RngStream rng[], int thread_id){
    return ((A) + ((B)-(A))*(rng[thread_id].RandU01()));
}

     /**-----**/
    /** Rod **/
   /**-----**/

/**
 Create a rod of a known length. The rod_no is an arbitrary identifier.
 Note that this creates a rod without initialising the contents of the
 arrays.
*/
Rod::Rod(int length, int set_rod_no):
    equil_r(new float[length]),
    equil_m(new float[length]),
    current_r(new float[length]),
    current_m(new float[length]),
    perturbed_x_energy_positive(new float[length]),
    perturbed_y_energy_positive(new float[length]),
    perturbed_z_energy_positive(new float[length]),
    twisted_energy_positive(new float[length]),
    perturbed_x_energy_negative(new float[length]),
    perturbed_y_energy_negative(new float[length]),
    perturbed_z_energy_negative(new float[length]),
    twisted_energy_negative(new float[length]),
    material_params(new float[length]),
    B_matrix(new float[length+(length/3)]),
    applied_forces(new float[length+(length/3)]),
    pinned_nodes(new bool[length/3])
    {}; 
        
/**
 Create a rod from a file. This won't do anything other than create the
 object - it won't even set up any arrays, because it doesn't know how
 long they'll have to be. After this, you have to call load_header and
 load_contents, which actually do the dirty work.
*/
Rod::Rod(std::string path, int set_rod_no):
    /** When we initialize from a file, we don't allocate arrays until we've loaded the file. **/
    line_start(0)
    {
        rod_no = set_rod_no;
    }; 
    
/**
 The contents of the rod, by default, are specified in SI units. Although
 it's possible to do everything in SI, you'll get more precision out of
 FFEA units. This function will convert all the units into FFEA units.
 When the file is written, the units are converted back automagically.
 The units are specified in mesoDimensions.h.
*/
Rod Rod::set_units(){
    /** Translate our units into the units specified in FFEA's mesoDimensions header file **/
    bending_response_factor = pow(mesoDimensions::length, 4)*mesoDimensions::pressure;
    spring_constant_factor = mesoDimensions::force/mesoDimensions::length;
    twist_constant_factor = mesoDimensions::force*mesoDimensions::length*mesoDimensions::length;
    
    /** And now the rod itself **/
    for (int i=0; i<length; i++){
        equil_r[i] /= mesoDimensions::length;
        equil_m[i] /= mesoDimensions::length;
        current_m[i] /= mesoDimensions::length;
        current_r[i] /= mesoDimensions::length;
        perturbed_x_energy_positive[i] /= mesoDimensions::Energy;
        perturbed_y_energy_positive[i] /= mesoDimensions::Energy;
        perturbed_z_energy_positive[i] /= mesoDimensions::Energy;
        twisted_energy_positive[i] /= mesoDimensions::Energy;
        perturbed_x_energy_negative[i] /= mesoDimensions::Energy;
        perturbed_y_energy_negative[i] /= mesoDimensions::Energy;
        perturbed_z_energy_negative[i] /= mesoDimensions::Energy;
        twisted_energy_negative[i] /= mesoDimensions::Energy;
        
        if (i%3 == 0){
            material_params[i] /= spring_constant_factor;
        }
        
        if (i%3 == 1){
            material_params[i] /= twist_constant_factor;
        }

        if (i%3 == 2){
            material_params[i] /= mesoDimensions::length;
        }
        
    }
    
    for (int i=0; i<length+length/3; i++){
        B_matrix[i] /= bending_response_factor;
    }
    
    return *this; /** Return a pointer to the object itself instead of void. Allows for method chaining! **/
}

     /**---------**/
    /** Updates **/
   /**---------**/

/**
 Do a timestep.
 This function contains two loops. Both are over the nodes. The first loop
 populates the contents of the energy arrays, which int we use to work out
 delta E. The second one uses those energies to compute dynamics and
 applies those dynamics to the position arrays.
*/
Rod Rod::do_timestep(RngStream rng[]){ // Most exciting method
    
    //The first loop is over all the nodes, and it computes all the energies for each one
    #pragma omp parallel for schedule(dynamic) //most of the execution time is spent in this first loop
    for (int node_no = 0; node_no<num_elements; node_no++){
        
        // if the node is pinned, we go to the next iteration of the loop (e.g. the next node)
        if (pinned_nodes[node_no] == true){
            continue;
        }

        // the cutoff values tell us how many nodes in our 5 node slice 'don't exist'
        // e.g. if we are at node i = n (at the end of the rod) the end cutoff will be 2
        int start_cutoff_val;
        int end_cutoff_val;
        int *start_cutoff = &start_cutoff_val;
        int *end_cutoff = &end_cutoff_val; // for the multiple return values
        set_cutoff_values(node_no, num_elements, start_cutoff, end_cutoff);
        
        // We need this e now because we need the previous value of e to do a material frame update
        // If you're curious about the [4][3] check out the get_perturbation_energy docs
        float p[4][3];
        load_p(p, current_r, node_no);
        
        float energies[3]; //bend, stretch, twist (temporary variable)
        
        // todo: these can have fewer arguments - maybe perturbation amount and cutoff values
        
        // We move the node backwards and forwards in each degree of freedom, so we end up calling get_perturbation_energy eight whole times
        // Fill the temporary variable with energies ( we basically pass the entire state of the rod to get_perturbation_energy)
        get_perturbation_energy( 
            perturbation_amount*0.5, //half one way, half the other
            x, // dimension (x, y, z are array indices, defined to be 1, 2, 3 at the top of this file, twist is = 4)
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies); // this is a void function with energies as the output

        // transfer them from the temporary variable to the real thing
        perturbed_x_energy_positive[node_no*3] = energies[stretch_index];
        perturbed_x_energy_positive[(node_no*3)+1] = energies[bend_index];
        perturbed_x_energy_positive[(node_no*3)+2] = energies[twist_index];   
                   
        get_perturbation_energy( //from rod_math
            perturbation_amount*0.5, 
            y, // dimension
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);
  
        perturbed_y_energy_positive[node_no*3] = energies[stretch_index];
        perturbed_y_energy_positive[(node_no*3)+1] = energies[bend_index];
        perturbed_y_energy_positive[(node_no*3)+2] = energies[twist_index];

        get_perturbation_energy( //from rod_math
            perturbation_amount*0.5,
            z, // dimension
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);
  
        perturbed_z_energy_positive[node_no*3] = energies[stretch_index];
        perturbed_z_energy_positive[(node_no*3)+1] = energies[bend_index];
        perturbed_z_energy_positive[(node_no*3)+2] = energies[twist_index];

        float twist_perturbation = 0.006283185; // 2pi/1000
        get_perturbation_energy( //from rod_math
            twist_perturbation*0.5,
            4, // twist dimension = 4
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);

        twisted_energy_positive[node_no*3] = energies[stretch_index];
        twisted_energy_positive[(node_no*3)+1] = energies[bend_index];
        twisted_energy_positive[(node_no*3)+2] = energies[twist_index];
        
        get_perturbation_energy( 
            perturbation_amount*-0.5,
            x, // dimension
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);

        perturbed_x_energy_negative[node_no*3] = energies[stretch_index];
        perturbed_x_energy_negative[(node_no*3)+1] = energies[bend_index];
        perturbed_x_energy_negative[(node_no*3)+2] = energies[twist_index];   
                   
        get_perturbation_energy( //from rod_math
            perturbation_amount*-0.5,
            y, // dimension
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);
            
        perturbed_y_energy_negative[node_no*3] = energies[stretch_index];
        perturbed_y_energy_negative[(node_no*3)+1] = energies[bend_index];
        perturbed_y_energy_negative[(node_no*3)+2] = energies[twist_index];

        get_perturbation_energy( //from rod_math
            perturbation_amount*-0.5,
            z, // dimension
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);
        
        perturbed_z_energy_negative[node_no*3] = energies[stretch_index];
        perturbed_z_energy_negative[(node_no*3)+1] = energies[bend_index];
        perturbed_z_energy_negative[(node_no*3)+2] = energies[twist_index];

        get_perturbation_energy( //from rod_math
            twist_perturbation*-0.5,
            4, // twist dimension = 4
            B_matrix,
            material_params,
            start_cutoff_val,
            end_cutoff_val,
            node_no,
            current_r,
            equil_r,
            current_m,
            equil_m,
            energies);
            
        twisted_energy_negative[node_no*3] = energies[stretch_index];
        twisted_energy_negative[(node_no*3)+1] = energies[bend_index];
        twisted_energy_negative[(node_no*3)+2] = energies[twist_index];
        
    }
    
    //This loop is for the dynamics
    for (int node_no = 0; node_no<num_elements; node_no++){
        
        // If the node is pinned, there's nothing to do
        if (pinned_nodes[node_no] == true){
            continue;
        }
        
        // Grab thread ID from openMP (needed for RNG)
        #ifdef USE_OPENMP
            int thread_id = omp_get_thread_num();
        #else
            int thread_id = 0;
        #endif
        
        // Get friction, needed for delta r and delta theta
        float translational_friction = get_translational_friction(this->viscosity, material_params[(node_no*3)+2], false);
        float length_for_friction = (get_absolute_length_from_array(equil_r, node_no, this->length) + get_absolute_length_from_array(equil_r, node_no-1, this->length))/2;
        float rotational_friction = get_rotational_friction(this->viscosity, material_params[(node_no*3)+2], length_for_friction, true);
        
        // Need these again
        float twist_perturbation = 0.006283185; // 2pi/1000

        // The material frame update requires that we grab the ith segment as it was before we did any dynamics           
        float previous_p_i[3];
        float r_i[3];
        float r_ip1[3];
        r_i[0] = current_r[node_no*3]; r_i[1] = current_r[(node_no*3)+1];  r_i[2] = current_r[(node_no*3)+2]; 
        r_ip1[0] = current_r[(node_no*3)+3]; r_ip1[1] = current_r[(node_no*3)+4];  r_ip1[2] = current_r[(node_no*3)+5];
        get_p_i(r_i, r_ip1, previous_p_i);
    
        // Get fluctuating force
        float x_noise = get_noise(timestep, kT, translational_friction, random_number(-0.5, 0.5, rng, thread_id));
        float y_noise = get_noise(timestep, kT, translational_friction, random_number(-0.5, 0.5, rng, thread_id));
        float z_noise = get_noise(timestep, kT, translational_friction, random_number(-0.5, 0.5, rng, thread_id));
        float twist_noise = get_noise(timestep, kT, rotational_friction, random_number(-0.5, 0.5, rng, thread_id));            

        // Sum our energies and use them to compute the force            
        float x_force = (perturbed_x_energy_negative[node_no*3]+perturbed_x_energy_negative[(node_no*3)+1]+perturbed_x_energy_negative[(node_no*3)+2] - (perturbed_x_energy_positive[node_no*3]+perturbed_x_energy_positive[(node_no*3)+1]+perturbed_x_energy_positive[(node_no*3)+2] ))/perturbation_amount;
        float y_force = (perturbed_y_energy_negative[node_no*3]+perturbed_y_energy_negative[(node_no*3)+1]+perturbed_y_energy_negative[(node_no*3)+2] - (perturbed_y_energy_positive[node_no*3]+perturbed_y_energy_positive[(node_no*3)+1]+perturbed_y_energy_positive[(node_no*3)+2] ))/perturbation_amount;
        float z_force = (perturbed_z_energy_negative[node_no*3]+perturbed_z_energy_negative[(node_no*3)+1]+perturbed_z_energy_negative[(node_no*3)+2] - (perturbed_z_energy_positive[node_no*3]+perturbed_z_energy_positive[(node_no*3)+1]+perturbed_z_energy_positive[(node_no*3)+2] ))/perturbation_amount;
        float twist_force = (twisted_energy_negative[node_no*3]+twisted_energy_negative[(node_no*3)+1]+twisted_energy_negative[(node_no*3)+2] - (twisted_energy_positive[node_no*3]+twisted_energy_positive[(node_no*3)+1]+twisted_energy_positive[(node_no*3)+2] ))/twist_perturbation;
                
        // Get applied force, if any
        float applied_force_x = applied_forces[node_no*4];
        float applied_force_y = applied_forces[(node_no*4)+1];
        float applied_force_z = applied_forces[(node_no*4)+2];
        float applied_force_twist = applied_forces[(node_no*4)+3];
        
        // Get delta r and delta twist
        float delta_r_x = get_delta_r(translational_friction, timestep, x_force, x_noise, applied_force_x); //from rod_math
        float delta_r_y = get_delta_r(translational_friction, timestep, y_force, y_noise, applied_force_y);
        float delta_r_z = get_delta_r(translational_friction, timestep, z_force, z_noise, applied_force_z);
        float delta_twist = get_delta_r(rotational_friction, timestep, twist_force, twist_noise, applied_force_twist);

        // Apply our delta x
        current_r[node_no*3] += delta_r_x;
        current_r[(node_no*3)+1] += delta_r_y;
        current_r[(node_no*3)+2] += delta_r_z;
        
        // A wee sanity check to stop your simulations from exploding horribly
        for (int i=0; i<length; i++){
            if (current_r[i] >= 30000000){
                std::cout << "node " << node_no << " frame " << frame_no << "\n";
                std::cout << "dynamics: " << delta_r_x << ", " << delta_r_y << ", " << delta_r_z << "\n";
                assert(current_r[i] < 30000000 && "r went crazy after initializing.");
            }
        }
        
        // If we're applying delta twist, we must load our new p_i back in
        if (node_no != num_elements - 1){ // The last node has no p_i, so it can't rotate
            float m_to_rotate[3];
            m_to_rotate[0] = current_m[node_no*3]; m_to_rotate[1] = current_m[(node_no*3)+1]; m_to_rotate[2] = current_m[(node_no*3)+2]; // take the relevant info out of the data structure
            float p_i[3];
            float r_i[3];
            float r_ip1[3]; // we need all of these from the rod again
            r_i[0] = current_r[node_no*3]; r_i[1] = current_r[(node_no*3)+1]; r_i[2] = current_r[(node_no*3)+2];
            r_ip1[0] = current_r[(node_no*3)+3]; r_ip1[1] = current_r[(node_no*3)+1+3]; r_ip1[2] = current_r[(node_no*3)+2+3];
            get_p_i(r_i, r_ip1, p_i); //from rod_math
            rodrigues_rotation(m_to_rotate, p_i, delta_twist, m_to_rotate); // work out the actual rotated value
            current_m[node_no*3] = m_to_rotate[0]; current_m[(node_no*3)+1] = m_to_rotate[1]; current_m[(node_no*3)+2] = m_to_rotate[2]; //put it back in the data structure
        }
        
        // If the element has moved, we need to update the material frame to have moved accordingly
        if (node_no != num_elements - 1){ // for last node index, material frame doesn't exist!
            float current_p_i[3]; 
            float m_to_fix[3];
            float m_i_prime[3];
            float r_i[3];
            float r_ip1[3]; //now: grab the quantities we need out of the data structure
            r_i[0] = current_r[node_no*3]; r_i[1] = current_r[(node_no*3)+1]; r_i[2] = current_r[(node_no*3)+2];
            r_ip1[0] = current_r[(node_no*3)+3]; r_ip1[1] = current_r[(node_no*3)+1+3]; r_ip1[2] = current_r[(node_no*3)+2+3];
            for (int i=0; i<3; i++){
                current_p_i[i] = r_ip1[i] - r_i[i];
            }
            m_to_fix[0] = current_m[node_no*3]; m_to_fix[1] = current_m[(node_no*3)+1]; m_to_fix[2] = current_m[(node_no*3)+2];
            update_m1_matrix(m_to_fix, previous_p_i, current_p_i, m_i_prime); //from rod_math - notice we're using the previous p_i we grabbed at the start of the function
            current_m[node_no*3] = m_i_prime[0]; current_m[(node_no*3)+1] = m_i_prime[1]; current_m[(node_no*3)+2] = m_i_prime[2]; // back into the data structure you go
        }

        step_no += 1; //we just did one timestep so increment this

    }
    
    return *this; 
}
    
     /**----**/
    /** IO **/
   /**----**/
    
/**
 Load the header info from a .rodtraj file (that's everything before the
 ---END HEADER--- line). Not all the info is read, some of it is for
 clarity. This populates some rod variables:
 - length - total length of the array (normally 3x the number of nodes)
 - num_elements - number of nodes in the rod
 - num_rods - number of rods in the simulation. Not used right now.
 - line_start - number of the line at which the trajectory begins. This
   variable is used by load_contents later on, to skip the header.
 - version - which version of the algorithm is this made for? 
 This method will also allocate the memory for all the arrays in the
 rod. Descriptions of those array are in rod_structure.h. Finally, it
 sets some default values for global simulation parameters. Eventually,
 these will be overwritten by parameters from the .ffea file.
*/
Rod Rod::load_header(std::string filename){
    rod_filename = filename;
    file_ptr = fopen(filename.c_str(),"a");
    
    /** This string denotes where the header info ends */
    const std::string rod_connections = "CONNECTIONS,ROD,0";
    
    /** Check that we can load our input file */
    std::ifstream infile(filename);
    if(!infile){std::cout << "Walrus IO error. Does input file exist?" << std::endl;}
    
    /** Iterate through file up to rod connections marker. */
    int n = 0;
    bool length_set = false; /** Use this to prevent rods with bad header info from being initialized **/
    for( std::string line; getline( infile, line ); ){
        if (n == 0){assert(line == "format,ffea_rod");} /** Check that format is valid FFEA_rod */
        if (n > 0 and line != rod_connections){
            
            /** Extract data from lines and deposit it into object */
            std::vector<std::string> line_vec;
            boost::split(line_vec, line, boost::is_any_of(","));
            
            /** Read in the contents */
            if (line_vec[0] == "version"){ this->rod_version = std::stod(line_vec[1]);}
            if (line_vec[0] == "length"){ this->length = std::stoi(line_vec[1]); length_set=true; }
            if (line_vec[0] == "num_elements"){ this->num_elements = std::stoi(line_vec[1]); }
            if (line_vec[0] == "num_rods"){ this->num_rods = std::stoi(line_vec[1]); }
        }   
        if (line == rod_connections){
            this->line_start = n; /** Set this variable so other methods know we've read the header */
            }
        n++;
    }
    
    assert(length_set == true && "Length of the rod has not been set.");
    
    /** Warn the user if there is a file version mismatch */
    if (fabs(this->rod_version - rod_software_version) > 0.0000001){
        std::cout << "Serious warning: The rod input was generated by a different version of the software to the one you are using. \n";
        std::cout << "Rod file version: " << this->rod_version << ", software version: " << rod_software_version << ".\n";
    }
    
    /** Now that we know the rod length, we can allocate the memory **/
    equil_r = static_cast<float *>(malloc(sizeof(float) * length));
    equil_m = static_cast<float *>(malloc(sizeof(float) * length));
    current_r = static_cast<float *>(malloc(sizeof(float) * length));
    current_m = static_cast<float *>(malloc(sizeof(float) * length));
    perturbed_x_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
    perturbed_y_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
    perturbed_z_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
    twisted_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
    perturbed_x_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
    perturbed_y_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
    perturbed_z_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
    twisted_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
    material_params = static_cast<float *>(malloc(sizeof(float) * length));
    B_matrix = static_cast<float *>(malloc(sizeof(float) * (length+(length/3)) ));
    applied_forces = static_cast<float *>(malloc(sizeof(float) * (length+(length/3)) ));
    pinned_nodes = static_cast<bool *>(malloc(sizeof(bool) * length/3));
    
    for (int i=0; i<length/3; i++){
        pinned_nodes[i] = false;
    }
    
    for (int i=0; i<length+(length/3); i++){
        applied_forces[i] = 0;
    }
    
    mersenne_twister = new std::mt19937(3535.442464);
    real_distribution = new std::uniform_real_distribution<float>(-0.5, 0.5);
    
    // These are hardcoded default values (temporary?) - eventualy they will load in from the .ffea script!
    viscosity_constant_factor = mesoDimensions::pressure*mesoDimensions::time; ///poiseuille
    this->viscosity = 0.6913*pow(10, -3)/viscosity_constant_factor;
    this->timestep = 1e-12/mesoDimensions::time;
    this->kT = 0;
    this->perturbation_amount = 0.001*pow(10,-9)/mesoDimensions::length; // todo: set this dynamically, maybe 1/1000 equilibrium length?
    
    return *this;
}

/**
 Add a constant force that acts upon the rod every timestep.
 Parameters: force[4], a 4-element array that specifies a 3-D force
 vector, with the last element being a torque instead.
 Node_index: the index of the node to apply these forces on.
 Note: to remove the force, you must call this function again with 0s,
 or it will continue appyling the force.
*/
Rod Rod::add_force(float force[4], int node_index){
    this->applied_forces[node_index*4] = force[0]/mesoDimensions::force;
    this->applied_forces[(node_index*4)+1] = force[1]/mesoDimensions::force;
    this->applied_forces[(node_index*4)+2] = force[2]/mesoDimensions::force;
    this->applied_forces[(node_index*4)+3] = force[3]/(mesoDimensions::force*mesoDimensions::length);
    return *this;
}
    
/**
 Load the current state of the rod. This function expects that load_header
 has already been called. This populates all of the already-initialised
 arrays containing the state of the rod. Note that it only contains the
 current state of the rod - the FFEA_rod python class is the only one
 that loads the rod trajectory.
*/
Rod Rod::load_contents(std::string filename){
    
    /** Make sure this method isn't called before loading header info */
    assert(line_start != 0 && "Rod header\rod file not found."); 
    
    std::ifstream infile(filename);
    int n = 0;
    int line_of_last_frame = 0;
    
    /** Get the line number of the last frame */
    for( std::string line; getline( infile, line ); ){
        std::vector<std::string> line_vec;
        boost::split(line_vec, line, boost::is_any_of(" "));
        if (line_vec[0] == "FRAME")
            { line_of_last_frame = n; }
        n++;
    }
    
    /** Check that we got the last frame */
    assert(line_of_last_frame != 0);
    
    /** Seek back to the start of the file */
    infile.clear();
    infile.seekg(0, std::ios::beg);
    n = 0;
    
    for( std::string line; getline( infile, line ); ){
        /** Find the last frame from the line number we got earlier */
        if (n > line_of_last_frame){
            /** Convert each string into a vector<string> */
            std::vector<std::string> line_vec;
            boost::split(line_vec, line, boost::is_any_of(","));
            /** Then convert that into a vector<float> */
            std::vector<float> line_vec_float;
            line_vec_float = stof_vec(line_vec);
            
            /** Check we're not going to overflow and ruin someone's life when we write into the array*/
            assert( (unsigned)length == line_vec_float.size() || (unsigned)length+(length/3) == line_vec_float.size() ); //it's definitely fine to cast length
            
            /** Set our rod data arrays to the raw .data() from the vector. */
            if (n == line_of_last_frame+1){ for (int i=0; i<length; i++) equil_r[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+2){ for (int i=0; i<length; i++) equil_m[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+3){ for (int i=0; i<length; i++) current_r[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+4){ for (int i=0; i<length; i++) current_m[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+5){ for (int i=0; i<length; i++) perturbed_x_energy_positive[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+6){ for (int i=0; i<length; i++) perturbed_y_energy_positive[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+7){ for (int i=0; i<length; i++) perturbed_z_energy_positive[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+8){ for (int i=0; i<length; i++) twisted_energy_positive[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+9){ for (int i=0; i<length; i++) perturbed_x_energy_negative[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+10){ for (int i=0; i<length; i++) perturbed_y_energy_negative[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+11){ for (int i=0; i<length; i++) perturbed_z_energy_negative[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+12){ for (int i=0; i<length; i++) twisted_energy_negative[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+13){ for (int i=0; i<length; i++) material_params[i] = line_vec_float.data()[i];}
            if (n == line_of_last_frame+14){ for (int i=0; i<length+(length/3); i++) B_matrix[i] = line_vec_float.data()[i];}

        }
        n++;
    }
    return *this;
}
    
/**
 Write the current state of the rod to a file specified by the pointer
 *file_ptr. This will convert from MesoDimensions to SI units, so if your
 values are already in SI units, they'll be wrong. 
*/
Rod Rod::write_frame_to_file(){
    this->frame_no += 1;
    std::fprintf(file_ptr, "FRAME %i ROD %i\n", frame_no, rod_no);
    write_array(equil_r, length, mesoDimensions::length);
    write_array(equil_m, length, mesoDimensions::length);
    write_array(current_r, length, mesoDimensions::length);
    write_array(current_m, length, mesoDimensions::length);
    write_array(perturbed_x_energy_positive, length, mesoDimensions::Energy);
    write_array(perturbed_y_energy_positive, length, mesoDimensions::Energy);
    write_array(perturbed_z_energy_positive, length, mesoDimensions::Energy);
    write_array(twisted_energy_positive, length, mesoDimensions::Energy);
    write_array(perturbed_x_energy_negative, length, mesoDimensions::Energy);
    write_array(perturbed_y_energy_negative, length, mesoDimensions::Energy);
    write_array(perturbed_z_energy_negative, length, mesoDimensions::Energy);
    write_array(twisted_energy_negative, length, mesoDimensions::Energy);
    write_mat_params_array(material_params, length, spring_constant_factor, twist_constant_factor, mesoDimensions::length);
    write_array(B_matrix, length+(length/3), bending_response_factor );
    return *this;
}
    
/**
 Write a single array to a file in the CSV format.
 Parameters:
  - *array_ptr - the pointer to the array that is to be written.
  - array_len - the length of the array.
  - unit_scale_factor - the unit conversion from the internal FFEA units
    to SI units.
*/
Rod Rod::write_array(float *array_ptr, int array_len, float unit_scale_factor){
    for (int i=0; i<array_len; i++){
        if (i<array_len-1){
            std::fprintf(file_ptr, "%e,", array_ptr[i]*unit_scale_factor);
        }
        else{
            std::fprintf(file_ptr, "%e", array_ptr[i]*unit_scale_factor);
        }
    }
    std::fprintf(file_ptr, "\n");
    return *this;
}

/**
 This function is almost identical to the one above, but it appllies
 different scale factors for objects in the array,
*/
Rod Rod::write_mat_params_array(float *array_ptr, int array_len, float stretch_scale_factor, float twist_scale_factor, float length_scale_factor){
    float scale_factors[3] = {stretch_scale_factor, twist_scale_factor, length_scale_factor};
    for (int i=0; i<array_len; i++){
        if (i<array_len-1){
            std::fprintf(file_ptr, "%e,", array_ptr[i]*scale_factors[i%3]);
        }
        else{
            std::fprintf(file_ptr, "%e", array_ptr[i]*scale_factors[i%3]);
        }
    }
    std::fprintf(file_ptr, "\n");
    return *this;
}

/**
 Close the previous file and create a new file, assigning that to the rod
 variable *file_ptr. This will also copy the contents of the previous
 file into this one.
*/
Rod Rod::change_filename(std::string new_filename){
    /** Get contents of current file */
    std::string current_file_contents = "";
    std::ifstream infile(this->rod_filename);
    for(int i=0 ; infile.eof()!=true ; i++){
        current_file_contents += infile.get();
    }
    current_file_contents.erase(current_file_contents.end()-1);
    infile.close();
    
    /** Create new file **/
    std::ofstream outfile(new_filename);
    outfile << current_file_contents;
    outfile.close();
    
    /** Update member variables (filename, pointer etc) */
    this->rod_filename = new_filename;
    fclose(file_ptr);
    file_ptr = fopen(new_filename.c_str(),"a");
    return *this;
}

/**
 Run the simulation for an arbitrary amount of time. If you start a
 rod exactly in its equilibrium state, chances are it's not going to be
 equilibrated, which can throw off some tests. It runs for a totally
 arbitrary 1e-7 seconds and does not save the trajectory from the
 equilibration.
*/
Rod Rod::equilibrate_rod(RngStream rng[]){
    int no_steps = 1e-7/timestep; // this is arbitrary
    for (int i=0; i<no_steps; i++){
        this->do_timestep(rng);
    }
    return *this;
}

/**
 Translate every node in the rod by a given translation vector, 
 translation_vec. The parameter float* r is the pointer to any array of
 node positions, e.g. this->current_r or this->equil_r. No return values,
 it just updates those arrays
 * */
Rod Rod::translate_rod(float* r, float translation_vec[3]){
    for(int i=0; i<this->length; i+=3){
        r[i] += translation_vec[0];
        r[i+1] += translation_vec[1];
        r[i+2] += translation_vec[2];
    }
    return *this;
}

/**
 Rotates the rod by the euler angles alpha, beta and gamma (or x, y
 and z if you prefer. This will update all the node positions AND the
 material frames will rotate as well. The rotations happen relative to
 each centroid, so if current_r and equil_r have different centroids,
 they will be rotated about different points.
*/
Rod Rod::rotate_rod(float euler_angles[3]){
    /** Put rod centroid on 0,0,0 */
    float equil_centroid[3];
    get_centroid(this->equil_r, this->length, equil_centroid);
    float equil_to_translate[3] = {-equil_centroid[0], -equil_centroid[1], -equil_centroid[2]};
    this->translate_rod(this->equil_r, equil_to_translate);
    
    float current_centroid[3];
    get_centroid(this->current_r, this->length, current_centroid);
    float current_to_translate[3] = {-current_centroid[0], -current_centroid[1], -current_centroid[2]};
    this->translate_rod(this->current_r, current_to_translate);
    
    /** Construct rotation matrix from euler angles **/
    float xe = euler_angles[0];
    float ye = euler_angles[1];
    float ze = euler_angles[2];
    
    float Rx[9] = {1, 0, 0, 0, cosf(xe), -sinf(xe), 0, sinf(xe), cosf(xe)};
    float Ry[9] = {cosf(ye), 0, sinf(ye), 0, 1, 0, -sinf(ye), 0, cosf(ye)};
    float Rz[9] = {cosf(ze), -sinf(ze), 0, sinf(ze), cosf(ze), 0, 0, 0, 1};
    
    float RyRx[9];
    matmul_3x3_3x3(Ry, Rx, RyRx);
    float RzRyRx[9];
    matmul_3x3_3x3(Rz, RyRx, RzRyRx);
        
    /** Apply rotation matrix **/
    for(int i=0; i<length; i+=3){
        float temp_vec[3];
        /** do not try to improve this */
        apply_rotation_matrix(this->equil_r+i, RzRyRx, temp_vec); 
        equil_r[i] = temp_vec[0]; equil_r[i+1] = temp_vec[1]; equil_r[i+2] = temp_vec[2];
        apply_rotation_matrix(this->current_r+i, RzRyRx, temp_vec);
        current_r[i] = temp_vec[0]; current_r[i+1] = temp_vec[1]; current_r[i+2] = temp_vec[2];
        apply_rotation_matrix(this->equil_m+i, RzRyRx, temp_vec);
        equil_m[i] = temp_vec[0]; equil_m[i+1] = temp_vec[1]; equil_m[i+2] = temp_vec[2];
        apply_rotation_matrix(this->current_m+i, RzRyRx, temp_vec);
        current_m[i] = temp_vec[0]; current_m[i+1] = temp_vec[1]; current_m[i+2] = temp_vec[2];
    }
    
    /** Move centroids back */
    this->translate_rod(this->current_r, current_centroid);
    this->translate_rod(this->equil_r, equil_centroid);
    
    return *this;
}

} //end namespace
