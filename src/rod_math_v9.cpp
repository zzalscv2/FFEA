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
 *      rod_math_v9.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk
 */

#define OUT ///< This is used to denote when a function modifies one of its parameters
#define _USE_MATH_DEFINES ///<  This has to come before including cmath

const static bool abort_on_fail = true; /// If the simulation becomes unstable, this will SIGABRT the whole program

#include <cmath>
#include <chrono>
#include <iostream>
#include <cstring>
#include <assert.h>
#include <random>
#include <stdexcept>
#include "dimensions.h"
#include <iomanip>    

namespace rod {
    
     /*---------*/
    /* Utility */
   /*---------*/

const double boltzmann_constant = 1.3806503e-23/mesoDimensions::Energy;

// A less weird way to access the contents of our arrays representing our vectors
static const int x = 0;
static const int y = 1;
static const int z = 2;

// For clarity, anytime you see an index [x], we're referring to the x
// dimension.

// Computing the energy for a perturbation requires four segments, defined
// by 5 nodes. They are i-2 to i+1, listed here.
static const int im2 = 0; ///< index of i-2nd thing
static const int im1 = 1; ///< index of i-1st thing
static const int i = 2; ///< index of ith thing
static const int ip1 = 3; ///< index of i+1th thing

#define OMP_SIMD_INTERNAL _Pragma("omp simd")

#define vec3d(x)for(int x = 0; x < 3; ++ x) ///< Shorthand to loop over elements of our 1d arrays representing 3d vectors

/** Hopefully I won't ever have to use this. */
void rod_abort(std::string message){
    std::cout << "There has been a cataclysmic error in FFEA_rod. Here it is:\n" << message << "\n";
    std::abort();
}

/**
 Check if a single value is simulation destroying. Here, simulation
 destroying means NaN or infinite.
*/
bool not_simulation_destroying(float x, std::string message){
    if(std::isnan(x) or std::isinf(x) ){
        if (abort_on_fail){ rod_abort(message); }
        else{ return false; }
    }
    return true;
}

/**
 This will do the same thing, but check an array 3 in length, and print
 a warning specifying which value it is.
*/
bool not_simulation_destroying(float x[3], std::string message){
    for (int i=0; i<3; i++){
        if (std::isnan(i) or std::isinf(i)){
            if (abort_on_fail){ rod_abort(message); }
            else{ return false; }
        }
    }
    return true;
}


/**
 Print the contents of an array to the stdout.
*/
void print_array(std::string array_name, float array[], int length){
    std::cout << array_name << " : [";
    for (int i = 0; i < length; i++){
        if (i != length-1){
            std::cout << array[i] << ", ";
        }
        else{
            std::cout << array[i];
        }
    }
    std::cout << "]\n";
}

// These are just generic vector functions that will be replaced by mat_vec_fns at some point

/**
 Normalize a 3-d vector. The there is no return value, but it populates
 an array whose pointer is specified as a function parameter, stl-style.
*/
void normalize(float in[3], OUT float out[3]){
    float absolute = sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]);
    vec3d(n){out[n] = in[n]/absolute;}
    not_simulation_destroying(out, "Noramlisation is simulation destroying."); // this is a cheaty way to give a c-style assertion an actual message.
}
/**
 Get the absolute value of a vector.
*/
float absolute(float in[3]){
    float absolute = sqrt(in[x]*in[x] + in[y]*in[y] + in[z]*in[z]);
    not_simulation_destroying(absolute, "Absolute value is simulation destroying.");
    assert(absolute>0 && "Absolute value is lower than zero (WHAT?).");
    return absolute;
}

/**
 Compute the cross product of a 3x1 vector x a 3x1 vector (the result is
 also a 3x1 vector).
*/
void cross_product(float a[3], float b[3], float out[3]){ // 3x1 x 3x1
    out[x] = (a[y]*b[z]) - (a[z] * b[y]);
    out[y] = (a[z]*b[x]) - (a[x] * b[z]);
    out[z] = (a[x]*b[y]) - (a[y] * b[x]);
    not_simulation_destroying(out, "Cross product is simulation destroying.");
}

/**
 Get the rotation matrix (3x3) that rotates a (3x1) onto b (3x1). 

  \f[R = I + [v]_\times + [v]_\times^2 \frac{1}{1+c}\f]
  Where
  \f[v = a \times b\f]
  \f[c = a \cdot b\f]
  \f[[v]_\times =   \begin{bmatrix}
    0 & -v_3 & v_2 \\
    v_3 & 0 & -v_1 \\
    -v_2 & v_1 & 0
    \end{bmatrix}\f]
 This seemed like the cheapest way to do it.
*/
void get_rotation_matrix(float a[3], float b[3], float rotation_matrix[9]){
    float v[3];
    cross_product(a,b,v);
    float c = (a[x]*b[x])+(a[y]*b[y])+(a[z]*b[z]);
    float vx[9];
    vx[0] = 0; vx[1] = -1*v[2]; vx[2] = v[1]; // vx = skew-symmetric cross product matrix
    vx[3] = v[2]; vx[4] = 0; vx[5] = -1*v[0];
    vx[6] = -1*v[1]; vx[7] = v[0]; vx[8] = 0;
    float m_f = 1/(1+c); // multiplication factor
    float identity_matrix[9] = {1,0,0,0,1,0,0,0,1};
    
    float vx_squared[9] = {-(v[1]*v[1])-(v[2]*v[2]), v[0]*v[1], v[0]*v[2], v[0]*v[1], -(v[0]*v[0])-(v[2]*v[2]), v[1]*v[2], v[0]*v[2], v[1]*v[2], -(-v[0]*-v[0])-(v[1]*v[1]) };
    
    for (int i=0; i<9; i++){
        rotation_matrix[i] = identity_matrix[i] + vx[i] + (vx_squared[i]*m_f);
    }
}

/**
 This is just a straight matrix multiplication, multiplyning the a column
 vector by a rotation matrix.
*/
void apply_rotation_matrix(float vec[3], float matrix[9], OUT float rotated_vec[3]){
    rotated_vec[0] = (vec[x]*matrix[0] + vec[y]*matrix[1] + vec[z]*matrix[2]);
    rotated_vec[1] = (vec[x]*matrix[3] + vec[y]*matrix[4] + vec[z]*matrix[5]);
    rotated_vec[2] = (vec[x]*matrix[6] + vec[y]*matrix[7] + vec[z]*matrix[8]);    
}

/**
 Same as above, but modifies a row vector instead of a row vector.
*/
void apply_rotation_matrix_row(float vec[3], float matrix[9], OUT float rotated_vec[3]){
    rotated_vec[0] = (vec[x]*matrix[0] + vec[y]*matrix[4] + vec[z]*matrix[7]);
    rotated_vec[1] = (vec[x]*matrix[2] + vec[y]*matrix[5] + vec[z]*matrix[8]);
    rotated_vec[2] = (vec[x]*matrix[3] + vec[y]*matrix[6] + vec[z]*matrix[9]);    
}

/**
 Dot product of two 3x3 matrices.
*/
void matmul_3x3_3x3(float a[9], float b[9], OUT float out[9]){
    out[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6]; out[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7]; out[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
    out[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6]; out[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7]; out[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
    out[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6]; out[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7]; out[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

// These are utility functions specific to the math for the rods

/**
 \f[ p_i = r_{i+1} - r_i \f]
  The segment \f$ e_i \f$ is the vector that runs from the node \f$ r_i \f$ to \f$ r_{i+1} \f$
*/
void get_p_i(float curr_r[3], float next_r[3], OUT float p_i[3]){
    vec3d(n){p_i[n] = next_r[n] - curr_r[n];}
    not_simulation_destroying(p_i, "Get_p_i is simulation destroying.");
}

/**
 \f[ {\mathbf  {v}}_{{\mathrm  {rot}}}={\mathbf  {v}}\cos \theta +({\mathbf  {k}}\times {\mathbf  {v}})\sin \theta +{\mathbf  {k}}({\mathbf  {k}}\cdot {\mathbf  {v}})(1-\cos \theta )~. \f]
 Where \f$ v_{rot} \f$ is the resultant vector, \f$ \theta \f$ is the angle to rotate,\f$ v \f$ is the original vector and \f$ k \f$ is the axis of rotation.
 This is Rodrigues' rotation formula, a cheap way to rotate a vector around an axis.
*/
void rodrigues_rotation(float v[3], float k[3], float theta, OUT float v_rot[3]){
    float k_norm[3]; 
    normalize(k, k_norm);
    float k_cross_v[3];
    float sin_theta = std::sin(theta); float cos_theta = std::cos(theta);
    cross_product(k_norm, v, k_cross_v);
    float right_multiplier = (1-cos_theta)*( (k_norm[x]*v[x])+(k_norm[y]*v[y])+(k_norm[z]*v[z]));
    float rhs[3];
    vec3d(n){ rhs[n] = right_multiplier * k_norm[n]; }
    vec3d(n){ v_rot[n] = cos_theta*v[n] + sin_theta*k_cross_v[n] + rhs[n]; }
    not_simulation_destroying(v_rot, "Rodrigues' rotation is simulation destroying.");
}

/**
 the c++ acos function will return nan for acos(>1), which we sometimes get (mostly 1.000001) due to
 some imprecisions. Obviously acos(a milllion) isn't a number, but for values very close to 1, we4
 will give the float the benefit of the doubt and say the acos is zero.
*/
float safe_cos(float in){
    if (in < 1.01){
        return 0;
    }
    
    float out = std::acos(in);
    if (std::isnan(out) or std::isinf(out)){
        if (abort_on_fail){ rod_abort("Cosine of something much larger than 1."); }
        else{ std::cout << "Warning: suspicious inverse cosine.\n"; return 0; }
    }
    
    return out;
}
 
 /**
 * Get the value of l_i.
 */
float get_l_i(float p_i[3], float p_im1[3]){
    return absolute(p_i) + absolute(p_im1);
}

     /*-----------------------*/
    /* Update Material Frame */
   /*-----------------------*/

/**
 \f[\widetilde{m_{1 i}}' = \widetilde{m_{1 i}} - ( \widetilde{m_{1 i}} \cdot \widetilde{l_i}) \widetilde{\hat{l_i}}\f]
 where \f$l\f$ is the normalized tangent, \f$m\f$ is the current material frame and \f$m'\f$ is the new one.
*/
void perpendicularize(float m_i[3], float p_i[3], OUT float m_i_prime[3]){
    float t_i[3];
    normalize(p_i, t_i);
    float m_i_dot_t_i = m_i[x]*t_i[x] + m_i[y]*t_i[y] + m_i[z]*t_i[z];
    vec3d(n){m_i_prime[n] = m_i[n] - m_i_dot_t_i*t_i[n] ;}
}

/**
 Say that the segment p_i is rotated into the position p_i_prime. This function rotates the material frame m_i
 by the same amount. Used to compute m_i of the 'perturbed' p_i values during the numerical differentation.
 And also when the new e_i values are computed at the end of each frame!
*/
void update_m1_matrix(float m_i[3], float p_i[3], float p_i_prime[3], float m_i_prime[3]){
    float rm[9];
    float m_i_rotated[3];
    get_rotation_matrix(p_i, p_i_prime, rm);
    apply_rotation_matrix(m_i, rm, m_i_rotated);
    perpendicularize(m_i_rotated, p_i, m_i_prime);
    normalize(m_i_prime, m_i_prime);
}

     /*------------------*/
    /* Compute Energies */
   /*------------------*/

/**
 \f[ E_{stretch} = \frac{1}{2}k(|\vec{p}_i| - |\widetilde{p}_i|)^2 \f]
 where \f$k\f$ is the spring constant, \f$p\f$ is the current segment and \f$m'\f$ is the equilbrium one.
*/
float get_stretch_energy(float k, float p_i[3], float p_i_equil[3]){
    
    float diff = absolute(p_i) - absolute(p_i_equil);
    
    float stretch_energy = (diff*diff*0.5*k)/absolute(p_i_equil);
    not_simulation_destroying(stretch_energy, "get_stretch_energy is simulation destroying.");
    
    return stretch_energy;
}

// todo: use OUT correctly on this fn

/**
 Use the previously defined rotation matrix functions to parallel transport a material frame
 m into the orientation m', from segment p_im1 to segment p_i.
*/
void parallel_transport(float m[3], float m_prime[3], float p_im1[3], float p_i[3]){
    float rm[9]; // rotation matrix
    get_rotation_matrix(p_im1, p_i, rm);
    apply_rotation_matrix(m, rm, m_prime);
}

/**
 \f[ E_{twist} = \frac{\beta}{l_i} \left( \Delta \theta_i - \Delta \widetilde{\theta}_i \right)^2 \f]
 Whereupon \f$l_i\f$ is \f$ |p_i| + |p_{i-1}| \f$, \f$\beta\f$ is the twisting energy constant, and
 \f[ \Delta\theta = \cos^{-1} ( P(m_{i+1}) \cdot m_i ) \f]
 Where P represents parallel transport.
*/
float get_twist_energy(float beta, float m_i[3], float m_im1[3], float m_i_equil[3], float m_im1_equil[3], float p_im1[3], float p_i[3], float p_im1_equil[3], float p_i_equil[3]){
    
    float l_i = get_l_i(p_im1_equil, p_i_equil);
    
    float p_i_norm[3];
    float p_im1_norm[3];
    float p_i_equil_norm[3];
    float p_im1_equil_norm[3];
    
    normalize(p_i, p_i_norm);
    normalize(p_im1, p_im1_norm);
    normalize(p_i_equil, p_i_equil_norm);
    normalize(p_im1_equil, p_im1_equil_norm);
    
    float m_prime[3];
    parallel_transport(m_im1, m_prime, p_im1_norm, p_i_norm);
    float m_equil_prime[3];
    parallel_transport(m_im1_equil, m_equil_prime, p_im1_equil_norm, p_i_equil_norm);
    float delta_theta = safe_cos( m_prime[x]*m_i[x] + m_prime[y]*m_i[y] + m_prime[z]*m_i[z] );
    float delta_theta_equil = safe_cos( m_equil_prime[x]*m_i_equil[x] + m_equil_prime[y]*m_i_equil[y] + m_equil_prime[z]*m_i_equil[z] );   
    float twist_energy = (beta/l_i)*((delta_theta - delta_theta_equil)*(delta_theta - delta_theta_equil));
    not_simulation_destroying(twist_energy, "get_twist_energy is simulation destroying.");

    return twist_energy;
}

/**
 \f[ \frac{2p_{i-1} \times p_i}{|p_i|\cdot|p_{i-1}| + p_{i-1}\cdot p_i } \f]
 Where \f$p_i\f$ and \f$p_{i-1}\f$ are the i-1 and ith segments, respectively.
*/
void get_kb_i(float p_im1[3], float p_i[3], OUT float kb_i[3]){
    float two_p_im1[3];
    vec3d(n){ two_p_im1[n] = p_im1[n] + p_im1[n];}
    float top[3];
    cross_product(two_p_im1, p_i, top);
    float bottom = (absolute(p_im1)*absolute(p_i))+( (p_im1[x]*p_i[x])+(p_im1[y]*p_i[y])+(p_im1[z]*p_i[z]));
    vec3d(n){ kb_i[n] = top[n]/bottom; }
    not_simulation_destroying(kb_i, "get_kb_i is simulation destroying.");
}

/**
 \f[ \omega(i,j) = \left( (k\vec{b})_i \cdot \vec{n}_j, -(k\vec{b})_i \cdot m_j \right)^T \f]
 Where \f$ (k\vec{b})_i \f$ is the curvature binormal, defined above, and \f$ m_j \f$ and \f$ n_j \f$ are the jth material axes.
*/
void get_omega_j_i(float kb_i[3], float n_j[3], float m_j[3], OUT float omega_j_i[2]){ //This is a column matrix not a vector
    omega_j_i[0] = (kb_i[x]*n_j[x])+(kb_i[y]*n_j[y])+(kb_i[z]*n_j[z]);
    omega_j_i[1] = -1*((kb_i[x]*m_j[x])+(kb_i[y]*m_j[y])+(kb_i[z]*m_j[z]));  
    not_simulation_destroying(omega_j_i[0], "get_omega_j_i is simulation destroying.");
    not_simulation_destroying(omega_j_i[1], "get_omega_j_i is simulation destroying.");
}

/**
 \f[ E_{bend} = \frac{1}{2 \widetilde{l}_i} \sum^i_{j=i-1} (\omega(i,j) - \widetilde{\omega}(i,j) )^T \widetilde{B}^i ( \omega(i,j) - \widetilde{\omega}(i,j) ) \f]
 Where \f$ \omega \f$ is the centreline curvature, defined above, \f$ B \f$ is the bending response matrix, and \f$l_i\f$ is \f$ |p_i| + |p_{i-1}| \f$
 
*/
float get_bend_energy(float omega_i_im1[2], float omega_i_im1_equil[2], float B_equil[4]){
    float delta_omega[2];
    delta_omega[0] = omega_i_im1[0] - omega_i_im1_equil[0];
    delta_omega[1] = omega_i_im1[1] - omega_i_im1_equil[1];
    float result = delta_omega[0]*(delta_omega[0]*B_equil[0] + delta_omega[1]*B_equil[2]) + delta_omega[1]*(delta_omega[0]*B_equil[1] + delta_omega[1]*B_equil[3]);
    not_simulation_destroying(result, "get_bend_energy is simulation destroying.");
    return result;
}

/**
 This function combines the curvature binormal, centerline curvature and bend energy formulae together, for a given set of segmments and material frames.
*/
float get_bend_energy_from_p(
        float p_im1[3],
        float p_i[3],
        float p_im1_equil[3],
        float p_i_equil[3],
        float n_im1[3],
        float m_im1[3],
        float n_im1_equil[3],
        float m_im1_equil[3],
        float n_i[3],
        float m_i[3],
        float n_i_equil[3],
        float m_i_equil[3],
        float B_i_equil[4],
        float B_im1_equil[4]){
            
    float p_i_norm[3];
    float p_im1_norm[3];
    float p_i_equil_norm[3];
    float p_im1_equil_norm[3];
    
    normalize(p_i, p_i_norm);
    normalize(p_im1, p_im1_norm);
    normalize(p_i_equil, p_i_equil_norm);
    normalize(p_im1_equil, p_im1_equil_norm);

    float l_i = get_l_i(p_i_equil, p_im1_equil);

    float kb_i[3];
    float kb_i_equil[3];
    get_kb_i(p_im1_norm, p_i_norm, kb_i);
    get_kb_i(p_im1_equil_norm, p_i_equil_norm, kb_i_equil);
    
    
    // Get omega and omega_equil for j = i-1
    float omega_j_im1[2];    
    get_omega_j_i(kb_i, n_im1, m_im1, omega_j_im1);
    
    float omega_j_im1_equil[2];    
    get_omega_j_i(kb_i_equil, n_im1_equil, m_im1_equil, omega_j_im1_equil);
        
    // And now for j = i
    float omega_j_i[2];    
    get_omega_j_i(kb_i, n_i, m_i, omega_j_i);
    
    float omega_j_i_equil[2];    
    get_omega_j_i(kb_i_equil, n_i_equil, m_i_equil, omega_j_i_equil);
    
    // Sum the bend energies between j = i-1 and j = i
    float bend_energy = 0;
    bend_energy += get_bend_energy(omega_j_i, omega_j_i_equil, B_i_equil);
    bend_energy += get_bend_energy(omega_j_im1, omega_j_im1_equil, B_im1_equil); //I THINK USING B_im1_equil IS WRONG
    bend_energy = bend_energy*(1/(2*l_i)); // constant!
    
    not_simulation_destroying(bend_energy, "get_bend_energy_from_p is simulation destroying.");
    
    if (bend_energy >= 1900000850){
        std::cout << "bend energy looks a bit large... here's a dump \n";
        print_array("p_im1", p_im1, 3);
        print_array("p_i", p_i, 3);
        print_array("p_im1_equil", p_im1_equil, 3);
        print_array("p_i_equil", p_i_equil, 3);
        print_array("n_im1_2", n_im1, 3);
        print_array("m_im1", m_im1, 3);
        print_array("n_im1_equil", n_im1_equil, 3);
        print_array("m_im1_equil", m_im1_equil, 3);
        
        print_array("n_i", n_i, 3);
        print_array("n_i_equil", n_i_equil, 3);
        print_array("m_i", m_i, 3);
        print_array("m_i_equil", m_i_equil, 3);
        print_array("B_i_equil", B_i_equil, 3);
        
        print_array("B_im1_equil", B_im1_equil, 4);
        std::cout << "l_equil = " << l_i << "\n";
        std::cout << "Energon Crystals = " << bend_energy << "\n";
//        assert(false);
    }
    
    return bend_energy;
}

void get_mutual_frame_inverse(float a[3], float b[3], OUT float mutual_frame[3]){
    float a_length = absolute(a);
    float b_length = absolute(b);
    vec3d(n){ mutual_frame[n] = (1/a_length)*a[n] + (1/b_length)*b[n]; }
    normalize(mutual_frame, mutual_frame);
}

float get_mutual_angle_inverse(float a[3], float b[3], float angle){
    float a_length = absolute(a);
    float b_length = absolute(b);
    float a_b_ratio = b_length/(b_length+a_length);
    return angle*a_b_ratio;
}

float get_bend_energy_mutual_parallel_transport(
        float p_im1[3],
        float p_i[3],
        float p_im1_equil[3],
        float p_i_equil[3],
        float n_im1[3],
        float m_im1[3],
        float n_im1_equil[3],
        float m_im1_equil[3],
        float n_i[3],
        float m_i[3],
        float n_i_equil[3],
        float m_i_equil[3],
        float B_i_equil[4],
        float B_im1_equil[4]){
    
    // get k_b
    float p_i_norm[3];
    float p_im1_norm[3];
    float p_i_equil_norm[3];
    float p_im1_equil_norm[3];
    
    normalize(p_i, p_i_norm);
    normalize(p_im1, p_im1_norm);
    normalize(p_i_equil, p_i_equil_norm);
    normalize(p_im1_equil, p_im1_equil_norm);

    float L_i = get_l_i(p_i_equil, p_im1_equil);

    float kb_i[3];
    float kb_i_equil[3];
    get_kb_i(p_im1_norm, p_i_norm, kb_i);
    get_kb_i(p_im1_equil_norm, p_i_equil_norm, kb_i_equil);
    
    // create our mutual l_i
    float mutual_l[3];
    float equil_mutual_l[3];
    get_mutual_frame_inverse(p_i, p_im1, OUT mutual_l);
    get_mutual_frame_inverse(p_i_equil, p_im1_equil, OUT equil_mutual_l);
    
    // parallel transport our existing material frames to our mutual l_i
    float m_im1_transported[3];
    float m_im1_equil_transported[3];
    parallel_transport(m_im1, m_im1_transported, p_im1_norm, mutual_l);
    parallel_transport(m_im1_equil, m_im1_equil_transported, p_im1_equil_norm, equil_mutual_l);
    
    float m_i_transported[3];
    float m_i_equil_transported[3];
    parallel_transport(m_i, m_i_transported, p_i_norm, mutual_l);
    parallel_transport(m_i_equil, m_i_equil_transported, p_i_equil_norm, equil_mutual_l);
    
    // get the angle between the two sets of material axes
    float angle_between_axes = safe_cos((m_i_transported[0] * m_im1_transported[0])+(m_i_transported[1] * m_im1_transported[1])+(m_i_transported[2] * m_im1_transported[2]));
    float angle_between_axes_equil = safe_cos((m_i_equil_transported[0] * m_im1_equil_transported[0])+(m_i_equil_transported[1] * m_im1_equil_transported[1])+(m_i_equil_transported[2] * m_im1_equil_transported[2]));

    float angle_to_rotate = get_mutual_angle_inverse(p_i, p_im1, angle_between_axes);
    float angle_to_rotate_equil = get_mutual_angle_inverse(p_i_equil, p_im1_equil, angle_between_axes_equil);
    
    float m_im1_rotated[3];
    float m_im1_rotated_equil[3];
    
    // rotate by the angle in question
    rodrigues_rotation(m_im1_transported, mutual_l, angle_to_rotate, m_im1_rotated);
    rodrigues_rotation(m_im1_equil_transported, equil_mutual_l, angle_to_rotate_equil, m_im1_rotated_equil);
    
    float n_im1_rotated[3];
    float n_im1_rotated_equil[3];
    
    cross_product(mutual_l, m_im1_rotated, n_im1_rotated);
    cross_product(equil_mutual_l, m_im1_rotated_equil, n_im1_rotated_equil);
    
    // finally get omega
    float omega_j_im1[2];    
    get_omega_j_i(kb_i, n_im1_rotated, m_im1_rotated, omega_j_im1);
    
    float omega_j_im1_equil[2];    
    get_omega_j_i(kb_i_equil, n_im1_rotated, m_im1_rotated_equil, omega_j_im1_equil);
    
    float bend_energy = 0;
    bend_energy += get_bend_energy(omega_j_im1, omega_j_im1_equil, B_i_equil);
    bend_energy = bend_energy*(1/(L_i)); // constant!

    not_simulation_destroying(bend_energy, "get_bend_energy_from_p is simulation destroying.");
    if(bend_energy >= 1900000850){ std::cout << "bend energy is very large. Please fire this up in gdb!\n"; }


    return bend_energy;
}

     /*----------*/
    /* Dynamics */
   /*----------*/
   
/**
 \f[ S_{translation} = \xi_{translation} = 6 \pi \mu a \f]
 Statement of Stokes law. The friction \f$ S \f$ for our dynamics can be computed from the viscosity of the medium, \f$\mu \f$, and the radius of the rod, \f$a\f$.
*/
float get_translational_friction(float viscosity, float radius, bool rotational){
    float friction = 6*M_PI*viscosity*radius;
    return friction;
}

/**
 \f[ S_{translation} = \xi_{translation} = 6 \pi \mu a \f]
 \f[ S_{rotation} = 8 \pi \mu a^3\f]
 \f[ S_{rotation} = 8 \pi \mu a^2 l \f]
 Both statements of Stokes law. The friction \f$ S \f$ for our dynamics can be computed from the viscosity of the medium, \f$\mu \f$, and the radius of the rod, \f$a\f$.
*/
float get_rotational_friction(float viscosity, float radius, float length, bool safe){
    if (safe){
        return 8*M_PI*viscosity*pow(radius, 2)*length;
    }
    else{
        return 8*M_PI*viscosity*pow(radius,3);
    }
}

/**
 \f[ F = \frac{\Delta E}{\Delta r} \f]
 Where \f$ F\f$ is force due to the energy gradient, \f$ \Delta E \f$ is the energy gradient, and \f$ \Delta r \f$ is the size of the perturbation.
 This is just a rearrangement of the work function.
*/
float get_force(float bend_energy, float stretch_energy, float delta_x){
    float result = bend_energy/delta_x + stretch_energy/delta_x;
    not_simulation_destroying(result, "get_force is simulation destroying.");
    return result;
}

/**
 \f[ T = \frac{\Delta E}{\Delta \theta} \f]
 Where \f$ T\f$ torque due to the energy gradient, \f$ \Delta E \f$ is the energy gradient, and \f$ \Delta \theta \f$ is the size of the perturbation.
 This is a rearrangement of the work function.
*/
float get_torque(float twist_energy, float delta_theta){
    float result = twist_energy/delta_theta;
    not_simulation_destroying(result, "get_torque is simulation destroying.");
    return result;
}

/**
 \f[ \Delta \underline{r}_i = \frac{\Delta t}{S} (F_c + F_{ext} + f_i) \f]
 Where \f$ \Delta \underline{r}_i  \f$ is the change in r (either x, y, z or \f$ \theta \f$, \f$ S \f$ is the viscous drag (derived from viscosity), \f$ F_C  \f$ is the force (or torque) due to the energy, \f$ F_{ext} \f$ is the external force being applied (if any) and \f$ f_i \f$ is the random force or torque.
 This expression is a rearrangement of the first order equation of motion for an object with viscous drag.
*/
float get_delta_r(float friction, float timestep, float force, float noise, float external_force){ // In a given dimension!
    float result = (timestep/friction)*(force + external_force + noise);
    not_simulation_destroying(result, "get_delta_x is simulation destroying.");
    return result;
}

/**
 \f[ g = \sqrt{ \frac{24k_B T \delta t }{S} } \f]
 Where \f$ g  \f$ is the force due to the thermal noise, \f$ T \f$ is the temperature of the system, and \f$ S \f$ is the viscous drag, derived from the viscosity.
 This expression is derived from the fluctuation dissipation relation for the langevin equation of the  .
*/
float get_noise(float timestep, float kT, float friction, float random_number){ // in a given dimension/DOF! 
    float result = std::sqrt( (24*kT*friction)/timestep );
    not_simulation_destroying(result*random_number, "get_noise is simulation destroying.");
    return result*random_number;
}

     /*------------*/
    /* Shorthands */
   /*------------*/

// Dear function pointer likers: I know, but I would rather directly call these functions.

/**
* This will load a region of a 1-d array containing nodes into an array
* of 4 p_i arrays, from i-2 to i+1. This is all the info we need to
* compute each type of energy.
*/   
void load_p(float p[4][3], float *r, int offset){
    int shift = (offset-2)*3;
    for (int j=0; j<4; j++){
        vec3d(n){ p[j][n] = r[shift+(j*3)+3+n] - r[shift+(j*3)+n]; }
    }
}

/**
 This does the same, only it loads m instead.
*/   
void load_m(float m_loaded[4][3], float *m, int offset){
    int shift = (offset-2)*3; // *3 for the 1-d array, -2 for offset 0 spanning i-2 to i+1
    for (int j=0; j<4; j++){
        m_loaded[j][0] = m[shift];
        m_loaded[j][1] = m[shift+1];
        m_loaded[j][2] = m[shift+2];
        shift += 3;
    }
}

/**
 This normalizes every segment in a 4-segment section of the rod.
*/   
void normalize_all(float p[4][3]){
    for (int j=0; j<4; j++){
        normalize(p[j], p[j]);
    }
}

/**
 This gets the absolute value of every segment in a 4-segment section of the rod.
*/   
void absolute_all(float p[4][3], float absolutes[4]){
    for (int j=0; j<4; j++){
        absolutes[j] = absolute(p[j]);
    }
}

/**
 This returns the value of m_2 (the cross product of e and m) for every
 segment in a 4-segment section of the rod.
*/   
void cross_all(float p[4][3], float m[4][3], OUT float n[4][3]){
    for (int j=0; j<4; j++){
        float p_norm[3];
        normalize(p[j], p_norm);
        cross_product(p_norm, m[j], n[j]);
    }
}

/**
 This computes the difference between two values of e for a given 4-segment
 section of the rod.
*/ 
void delta_e_all(float p[4][3], float new_p[4][3], OUT float delta_p[4][3]){
    for (int j=0; j<4; j++){
        delta_p[j][x] = new_p[j][x] - p[j][x];
        delta_p[j][y] = new_p[j][y] - p[j][y];
        delta_p[j][z] = new_p[j][z] - p[j][z]; 
    }
}

/**
 This performs the material frame update described earlier on every segment
 in a 4-segment section of the rod. Because this operation is more expensive
 than the others in this list, it uses a lookup table, and skips non-existent
 segments. 
*/ 
void update_m1_matrix_all(float m[4][3], float p[4][3], float p_prime[4][3], OUT float m_prime[4][3], int start_cutoff, int end_cutoff){
// I've tried writing 'clever' versions of this
// but ultimately it's clearer to just write the lookup table explicitly    
    if (start_cutoff == 0 and end_cutoff == 0){ //Somewhere in the middle of the rod
        update_m1_matrix(m[0], p[0], p_prime[0], m_prime[0]);
        update_m1_matrix(m[1], p[1], p_prime[1], m_prime[1]);
        update_m1_matrix(m[2], p[2], p_prime[2], m_prime[2]);
        update_m1_matrix(m[3], p[3], p_prime[3], m_prime[3]);
    }
    else if (start_cutoff == 1 and end_cutoff == 0){ // one node at the start is cut off
        update_m1_matrix(m[1], p[1], p_prime[1], m_prime[1]);
        update_m1_matrix(m[2], p[2], p_prime[2], m_prime[2]);
        update_m1_matrix(m[3], p[3], p_prime[3], m_prime[3]);
    }
    else if (start_cutoff == 2 and end_cutoff == 0){ // two at the start are cut off
        update_m1_matrix(m[2], p[2], p_prime[2], m_prime[2]);
        update_m1_matrix(m[3], p[3], p_prime[3], m_prime[3]);
    }
    else if (start_cutoff == 0 and end_cutoff == 1){ // one at the end is cut off
        update_m1_matrix(m[0], p[0], p_prime[0], m_prime[0]);
        update_m1_matrix(m[1], p[1], p_prime[1], m_prime[1]);
        update_m1_matrix(m[2], p[2], p_prime[2], m_prime[2]);
    }
    else if (start_cutoff == 0 and end_cutoff == 2){ // two at the end are cut off
        update_m1_matrix(m[0], p[0], p_prime[0], m_prime[0]);
        update_m1_matrix(m[1], p[1], p_prime[1], m_prime[1]);
    }
    else{
        throw std::invalid_argument("Length of the rod must be larger than 3");
    }
}

void load_B_all(float B[4][4], float *B_matrix, int offset){
    int shift = (offset-2)*4; // *3 for the 1-d array, -2 for offset 0 spanning i-2 to i+1
    for (int n=0; n<4; n++){
        B[n][0] = B_matrix[shift+0];
        B[n][1] = B_matrix[shift+1];
        B[n][2] = B_matrix[shift+2];
        B[n][3] = B_matrix[shift+3];
        shift += 4;
    }
}

     /*----------*/
    /* Utility  */
   /*----------*/

/**
 This is a utility function that will (arbitrarily) create a 4x4
 diagonal matrix with a symmetric bending response B.
*/ 
void make_diagonal_B_matrix(float B, OUT float B_matrix[4]){
    B_matrix[0] = B;
    B_matrix[1] = 0;
    B_matrix[2] = 0;
    B_matrix[3] = B;
}

// See notes on cutoff in get_perturbation_energy
/**
 The two variables, start and end cutoff, determine whether to skip some
 calculations (energies, material frame updates, etc). They're computed
 based on the position of the current node, and whether that node's
 4-segment 'window of influence' goes off the side of the rod or not.
 The values of start and end cutoff are how many nodes on their respective
 ends of the rod do not exist.
*/ 
void set_cutoff_values(int p_i_node_no, int num_nodes, OUT int *start_cutoff, int *end_cutoff){
    *start_cutoff = 0;
    *end_cutoff = 0;
    if (p_i_node_no == 0){*start_cutoff = 2;}
    if (p_i_node_no == 1){*start_cutoff = 1;}
    if (p_i_node_no == num_nodes-2){*end_cutoff = 1;}
    if (p_i_node_no == num_nodes-1){*end_cutoff = 2;}
}

/** Returns the absolute length of an element, given an array of elements
 *  and the index of the element itself.   */
float get_absolute_length_from_array(float* array, int node_no, int length){
    if (node_no*3 >= length){
        return 0;
    }
    else if (node_no*3 < 0){
        return 0;
    }
    else{
        float p[3] = {array[node_no*3], array[(node_no*3)+1], array[(node_no*3)+2]};
        return absolute(p);
    }
} 

/** Get the centroid of a particular rod, specified by the array of node
 positions for that rod (and the length). Updates the 'centroid' array
 given as a parameter. */
void get_centroid(float* r, int length, OUT float centroid[3]){
    float sum_pos[3] = {0,0,0};
    for(int i=0; i<length; i+=3){
        sum_pos[0] += r[i];
        sum_pos[1] += r[i+1];
        sum_pos[2] += r[i+2];
    }
    centroid[0] = sum_pos[0]/(length/3);
    centroid[1] = sum_pos[1]/(length/3);
    centroid[2] = sum_pos[2]/(length/3);
}

     /*-------------------------------*/
    /* Move the node, get the energy */
   /*-------------------------------*/

/**
 The get_perturbation_energy function ties together everything in this file.
 It will compute the energy in a specified degree of freedom for a given node.
   - perturbation_amount - the amount of perturbation to do in the numerical differentiation.
   - perturbation_dimension - which dimension to get dE/dr in (x,y,z or twist)
   - B_equil, k and beta - the bending response matrix, spring and twist constants
   - start and end cutoff, p_i_node_no - your position in the rod
   - r_all, r_all_equil, m_all, m_all_equil - pointers to arrays containing the complete state of the rod
   - energies - an array containing 3 values - stretch, bend and twist energy.

  It works as follows:
   - Load up a 2-d array for each e, m, e_equil and m_equil in the 4-segment zone needed to compute a dE\dr
   - Perturb a degree of freedom and update the material frame accordingly
   - Compute the value of j
   - Compute each energy based on the contribution from each segment affected by the numerical differentiation
*/ 
void get_perturbation_energy(
        float perturbation_amount,
        int perturbation_dimension,
        float *B_matrix,
        float *material_params,
        int start_cutoff,
        int end_cutoff,
        int p_i_node_no,
        float *r_all,
        float *r_all_equil,
        float *m_all,
        float *m_all_equil,
        OUT
        float energies[3]
    ){
        
    // Put a 5-node segment onto the stack.
    // We need to make a copy of it, because we'l be modifying it for our
    // Numerical differentiation later on.

    float B_equil[4][4];
    load_B_all(B_equil, B_matrix, p_i_node_no);
    
    // We'll end up modifying this, but we need the original later to update the material frame
    float original_p[4][3];   
    load_p(original_p, r_all, p_i_node_no);

    float p[4][3]; // the perturbed e
    float p_equil[4][3];
    float m[4][3];
    float m_equil[4][3];
    float material[4][3]; // 0 = k (stretch), 1 = beta (twist), 2 = unused (for now)
    
    // Compute e from m, and load into an appropriate data structure (a 2-D array)
    load_p(p, r_all, p_i_node_no);
    load_p(p_equil, r_all_equil, p_i_node_no);
    load_m(m, m_all, p_i_node_no);    
    load_m(m_equil, m_all_equil, p_i_node_no);
    load_m(material, material_params, p_i_node_no);
    
    // Apply our perturbation in x, y or z (for numerical differentiation)
    if (perturbation_dimension < 4 and perturbation_amount != 0){ //e.g. if we're perturbing x, y, or z
        p[im1][perturbation_dimension] += perturbation_amount;
        p[i][perturbation_dimension] -= perturbation_amount;
    }
    
    // If we perturb our angle instead, we apply a rodrigues rotation.
    if(perturbation_dimension == 4 and perturbation_amount != 0){ // if we're perturbing the twist
        rodrigues_rotation(m[i], p[i], perturbation_amount, m[i]);
    }

    // If we've perturbed it in x, y, or z, we need to update m, and then adjust it to make sure it's perpendicular
    if (perturbation_dimension < 4 and perturbation_amount != 0){ 
        update_m1_matrix_all(m, original_p, p, m, start_cutoff, end_cutoff);
    }

    // Normalize m1, just to be sure (the maths doesn't work if it's not normalized)
    normalize_all(m);
    normalize_all(m_equil);

    // Compute m_i_2 (we know it's perpendicular to e_i and m_i_1, so this shouldn't be too hard)
    float n[4][3];
    float n_equil[4][3];
    cross_all(p, m, n);
    cross_all(p_equil, m_equil, n_equil);
    
    // Compute unperturbed energy.
    // I could make this less verbose, but the explicit lookup table is a bit clearer about what's going on.
    // The basic idea is: if we're close to the 'edge' of the rod, don't compute energies for non-existent nodes! Because they are out of bounds!
    float bend_energy = 0;
    float stretch_energy = 0;
    float twist_energy = 0;
        
    if (start_cutoff == 0 and end_cutoff == 0){
        bend_energy += get_bend_energy_mutual_parallel_transport(p[im2], p[im1], p_equil[im2], p_equil[im1], n[im2], m[im2], n_equil[im2], m_equil[im2], n[im1], m[im1], n_equil[im1], m_equil[im1], B_equil[im1], B_equil[im2]);
        bend_energy += get_bend_energy_mutual_parallel_transport(p[im1], p[i], p_equil[im1], p_equil[i], n[im1], m[im1], n_equil[im1], m_equil[im1], n[i], m[i], n_equil[i], m_equil[i], B_equil[i], B_equil[im1]);
        stretch_energy += get_stretch_energy(material[im1][0], p[im1], p_equil[im1]);
        twist_energy += get_twist_energy(material[i][1], m[i], m[im1], m_equil[i], m_equil[im1], p[im1], p[i], p_equil[im1], p_equil[i]);
        stretch_energy += get_stretch_energy(material[i][0], p[i], p_equil[i]);
        twist_energy += get_twist_energy(material[ip1][1], m[ip1], m[i], m_equil[ip1], m_equil[i], p[i], p[ip1], p_equil[i], p_equil[ip1]);
        bend_energy += get_bend_energy_mutual_parallel_transport(p[i], p[ip1], p_equil[i], p_equil[ip1], n[i], m[i], n_equil[i], m_equil[i], n[ip1], m[ip1], n_equil[ip1], m_equil[ip1], B_equil[ip1], B_equil[i]);
        twist_energy += get_twist_energy(material[im1][1], m[im1], m[im2], m_equil[im1], m_equil[im2], p[im2], p[im1], p_equil[im2], p_equil[im1]);
    }
    else if (start_cutoff == 1 and end_cutoff == 0){
        bend_energy += get_bend_energy_mutual_parallel_transport(p[im1], p[i], p_equil[im1], p_equil[i], n[im1], m[im1], n_equil[im1], m_equil[im1], n[i], m[i], n_equil[i], m_equil[i], B_equil[i], B_equil[im1]);
        stretch_energy += get_stretch_energy(material[im1][0], p[im1], p_equil[im1]);
        twist_energy += get_twist_energy(material[i][1], m[i], m[im1], m_equil[i], m_equil[im1], p[im1], p[i], p_equil[im1], p_equil[i]);
        stretch_energy += get_stretch_energy(material[i][0], p[i], p_equil[i]);
        twist_energy += get_twist_energy(material[ip1][1], m[ip1], m[i], m_equil[ip1], m_equil[i], p[i], p[ip1], p_equil[i], p_equil[ip1]);
        bend_energy += get_bend_energy_mutual_parallel_transport(p[i], p[ip1], p_equil[i], p_equil[ip1], n[i], m[i], n_equil[i], m_equil[i], n[ip1], m[ip1], n_equil[ip1], m_equil[ip1], B_equil[ip1], B_equil[i]);
    }
    else if (start_cutoff == 2 and end_cutoff == 0){
        stretch_energy += get_stretch_energy(material[i][0], p[i], p_equil[i]);
        twist_energy += get_twist_energy(material[ip1][1], m[ip1], m[i], m_equil[ip1], m_equil[i], p[i], p[ip1], p_equil[i], p_equil[ip1]);
        bend_energy += get_bend_energy_mutual_parallel_transport(p[i], p[ip1], p_equil[i], p_equil[ip1], n[i], m[i], n_equil[i], m_equil[i], n[ip1], m[ip1], n_equil[ip1], m_equil[ip1], B_equil[ip1], B_equil[i]);
    }
    else if (start_cutoff == 0 and end_cutoff == 1){
        bend_energy += get_bend_energy_mutual_parallel_transport(p[im2], p[im1], p_equil[im2], p_equil[im1], n[im2], m[im2], n_equil[im2], m_equil[im2], n[im1], m[im1], n_equil[im1], m_equil[im1], B_equil[im1], B_equil[im2]);
        bend_energy += get_bend_energy_mutual_parallel_transport(p[im1], p[i], p_equil[im1], p_equil[i], n[im1], m[im1], n_equil[im1], m_equil[im1], n[i], m[i], n_equil[i], m_equil[i], B_equil[i], B_equil[im1]);
        stretch_energy += get_stretch_energy(material[im1][0], p[im1], p_equil[im1]);
        twist_energy += get_twist_energy(material[i][1], m[i], m[im1], m_equil[i], m_equil[im1], p[im1], p[i], p_equil[im1], p_equil[i]);
        stretch_energy += get_stretch_energy(material[i][0], p[i], p_equil[i]);
        twist_energy += get_twist_energy(material[im1][1], m[im1], m[im2], m_equil[im1], m_equil[im2], p[im2], p[im1], p_equil[im2], p_equil[im1]);
    }    
    else if (start_cutoff == 0 and end_cutoff == 2){
        bend_energy += get_bend_energy_mutual_parallel_transport(p[im2], p[im1], p_equil[im2], p_equil[im1], n[im2], m[im2], n_equil[im2], m_equil[im2], n[im1], m[im1], n_equil[im1], m_equil[im1], B_equil[im1], B_equil[im2]);
        stretch_energy += get_stretch_energy(material[im1][0], p[im1], p_equil[im1]);
        twist_energy += get_twist_energy(material[im1][1], m[im1], m[im2], m_equil[im1], m_equil[im2], p[im2], p[im1], p_equil[im2], p_equil[im1]);
     }
    else{
        throw std::invalid_argument("Length of the rod must be larger than 3");
    }

    energies[0] = bend_energy;
    energies[1] = stretch_energy;
    energies[2] = twist_energy;
}
//   _ _
//  (0v0)  I AM DEBUG OWL. PUT ME IN YOUR
//  (| |)  SOURCE CODE AND IT WILL BE BUG
//   W-W   FREE FOREVER. HOOT HOOT! 

} //end namespace
