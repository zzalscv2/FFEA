#include <math.h>
#include <stdlib.h>

static const int x = 0;
static const int y = 1;
static const int z = 2;
static const int im2 = 0; ///< index of i-2nd thing
static const int im1 = 1; ///< index of i-1st thing
static const int i = 2; ///< index of ith thing
static const int ip1 = 3; ///< index of i+1th thing
#define vec3d(x)for(int x = 0; x < 3; ++ x)
#define OUT
static const float rod_software_version = 0.3;

void normalize(float in[3], OUT float out[3]){
    float absolute = sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]);
    vec3d(n){out[n] = in[n]/absolute;}
}

void precise_normalize(float in[3], float out[3]){
    double in_double[3] = {(double)in[0], (double)in[1], (double)in[2]};
    float absolute = sqrt(in_double[0]*in_double[0] + in_double[1]*in_double[1] + in_double[2]*in_double[2]);
    float absolute_float = (float)absolute;
    vec3d(n){out[n] = in[n]/absolute_float;}

}

float absolute(float in[3]){
    float absolute = sqrt(in[x]*in[x] + in[y]*in[y] + in[z]*in[z]);
    return absolute;
}

void cross_product(float a[3], float b[3], float out[3]){ // 3x1 x 3x1
    out[x] = (a[y]*b[z]) - (a[z] * b[y]);
    out[y] = (a[z]*b[x]) - (a[x] * b[z]);
    out[z] = (a[x]*b[y]) - (a[y] * b[x]);
}

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

void apply_rotation_matrix(float vec[3], float matrix[9], OUT float rotated_vec[3]){
    rotated_vec[0] = (vec[x]*matrix[0] + vec[y]*matrix[1] + vec[z]*matrix[2]);
    rotated_vec[1] = (vec[x]*matrix[3] + vec[y]*matrix[4] + vec[z]*matrix[5]);
    rotated_vec[2] = (vec[x]*matrix[6] + vec[y]*matrix[7] + vec[z]*matrix[8]);    
}

float safe_cos(float in){
    
    float absin = abs(in);
    
    if (absin >= 1){
        return 0;
    }
    
    else{
        return acos(absin);
    }
}

float get_l_i(float p_i[3], float p_im1[3]){
    return (absolute(p_i) + absolute(p_im1))/2.0;
}
float get_signed_angle(float m1[3], float m2[3], float l[3]){
    float m2_cross_m1[3];
    cross_product(m2, m1, m2_cross_m1);
    return atan2( (m2_cross_m1[0] * l[0] + m2_cross_m1[1] * l[1] + m2_cross_m1[2] * l[2]), (m1[0]*m2[0]+m1[1]*m2[1]+m1[2]*m2[2]) );
}

float get_stretch_energy(float k, float p_i[3], float p_i_equil[3]){
    float diff = absolute(p_i) - absolute(p_i_equil);
    float stretch_energy = (diff*diff*0.5*k)/absolute(p_i_equil);
    return stretch_energy;
}
   

void parallel_transport(float m[3], float m_prime[3], float p_im1[3], float p_i[3]){
    float rm[9]; // rotation matrix
    get_rotation_matrix(p_im1, p_i, rm);
    apply_rotation_matrix(m, rm, m_prime);
}

float get_twist_energy(float beta, float m_i[3], float m_im1[3], float m_i_equil[3], float m_im1_equil[3], float p_im1[3], float p_i[3], float p_im1_equil[3], float p_i_equil[3]){
    
    float l_i = get_l_i(p_im1_equil, p_i_equil);
    
    float p_i_norm[3];
    float p_im1_norm[3];
    float p_i_equil_norm[3];
    float p_im1_equil_norm[3];
    
    float m_i_norm[3];
    float m_i_equil_norm[3];
    
    normalize(p_i, p_i_norm);
    normalize(p_im1, p_im1_norm);
    normalize(p_i_equil, p_i_equil_norm);
    normalize(p_im1_equil, p_im1_equil_norm);
    
    precise_normalize(m_i, m_i_norm);
    precise_normalize(m_i_equil, m_i_equil_norm);
    
    float m_prime[3];
    parallel_transport(m_im1, m_prime, p_im1_norm, p_i_norm);
    float m_equil_prime[3];
    parallel_transport(m_im1_equil, m_equil_prime, p_im1_equil_norm, p_i_equil_norm);
    
    precise_normalize(m_prime, m_prime);
    precise_normalize(m_equil_prime, m_equil_prime);
    
    float delta_theta = get_signed_angle(m_prime, m_i_norm, p_i_norm);
    float delta_theta_equil = get_signed_angle(m_equil_prime, m_i_equil_norm, p_i_equil_norm);
        
    float twist_energy = beta/(l_i*2) * pow(fmod( delta_theta - delta_theta_equil + M_PI, 2*M_PI) - M_PI, 2);

    return twist_energy;
}

void get_kb_i(float p_im1[3], float p_i[3], OUT float kb_i[3]){
    float two_p_im1[3];
    vec3d(n){ two_p_im1[n] = p_im1[n] + p_im1[n];}
    float top[3];
    cross_product(two_p_im1, p_i, top);
    float bottom = (absolute(p_im1)*absolute(p_i))+( (p_im1[x]*p_i[x])+(p_im1[y]*p_i[y])+(p_im1[z]*p_i[z]));
    vec3d(n){ kb_i[n] = top[n]/bottom; }
}

void get_omega_j_i(float kb_i[3], float n_j[3], float m_j[3], OUT float omega_j_i[2]){ //This is a column matrix not a vector
    omega_j_i[0] = (kb_i[x]*n_j[x])+(kb_i[y]*n_j[y])+(kb_i[z]*n_j[z]);
    omega_j_i[1] = -1*((kb_i[x]*m_j[x])+(kb_i[y]*m_j[y])+(kb_i[z]*m_j[z]));  
}

float get_bend_energy(float omega_i_im1[2], float omega_i_im1_equil[2], float B_equil[4]){
    float delta_omega[2];
    delta_omega[0] = omega_i_im1[0] - omega_i_im1_equil[0];
    delta_omega[1] = omega_i_im1[1] - omega_i_im1_equil[1];
    float result = delta_omega[0]*(delta_omega[0]*B_equil[0] + delta_omega[1]*B_equil[2]) + delta_omega[1]*(delta_omega[0]*B_equil[1] + delta_omega[1]*B_equil[3]);
    return result;
}

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
    
    
    float omega_j_im1[2];    
    get_omega_j_i(kb_i, n_im1, m_im1, omega_j_im1);
    
    float omega_j_im1_equil[2];    
    get_omega_j_i(kb_i_equil, n_im1_equil, m_im1_equil, omega_j_im1_equil);
        
    float omega_j_i[2];    
    get_omega_j_i(kb_i, n_i, m_i, omega_j_i);
    
    float omega_j_i_equil[2];    
    get_omega_j_i(kb_i_equil, n_i_equil, m_i_equil, omega_j_i_equil);
    
    float bend_energy = 0;
    bend_energy += get_bend_energy(omega_j_i, omega_j_i_equil, B_i_equil);
    bend_energy += get_bend_energy(omega_j_im1, omega_j_im1_equil, B_im1_equil); //I THINK USING B_im1_equil IS WRONG
    bend_energy = bend_energy*(1/(2*l_i)); // constant!
        
    return bend_energy;
}

float get_weights(float a[3], float b[3]){
    float a_length = absolute(a);
    float b_length = absolute(b);
    float weight1 = a_length/(a_length+b_length);
    return weight1; // weight2 = 1-weight1
}

void get_mutual_element_inverse(float pim1[3], float pi[3], float weight, OUT float mutual_element[3]){
    float pim1_norm[3];
    float pi_norm[3];
    normalize(pi, pi_norm);
    normalize(pim1, pim1_norm);
    vec3d(n){mutual_element[n] = (1/weight)*pim1_norm[n] + (1/(1-weight))*pi_norm[n]; }
    normalize(mutual_element, mutual_element);
}

void get_mutual_axes_inverse(float mim1[3], float mi[3], float weight, OUT float m_mutual[3]){
    float mi_length = absolute(mi);
    float mim1_length = absolute(mim1);
    vec3d(n){m_mutual[n] = (mim1[n]*(1.0/weight) + mi[n]*(1.0/(1-weight)))/(mi_length+mim1_length); }
    normalize(m_mutual, m_mutual);
}

float get_bend_energy_mutual_parallel_transport(
        float p_im1[3],
        float p_i[3],
        float p_im1_equil[3],
        float p_i_equil[3],
        float m_im1[3],
        float m_im1_equil[3],
        float m_i[3],
        float m_i_equil[3],
        float B_i_equil[4],
        float B_im1_equil[4],
        float omega[2],
        float omega_equil[2]
        ){

    // get k_b
    float p_i_norm[3];
    float p_im1_norm[3];
    float p_i_equil_norm[3];
    float p_im1_equil_norm[3];
    
    normalize(p_i, p_i_norm);
    normalize(p_im1, p_im1_norm);
    normalize(p_i_equil, p_i_equil_norm);
    normalize(p_im1_equil, p_im1_equil_norm);

    float n_im1[3];
    float n_i[3];
    float n_i_equil[3];
    float n_im1_equil[3];
    cross_product(m_im1, p_im1_norm, n_im1);
    cross_product(m_i, p_i_norm, n_i);
    cross_product(m_i_equil, p_i_equil_norm, n_i_equil);
    cross_product(m_im1_equil, p_im1_equil_norm, n_im1_equil);

    float L_i = get_l_i(p_i_equil, p_im1_equil);

    float kb_i[3];
    float kb_i_equil[3];
    get_kb_i(p_im1_norm, p_i_norm, kb_i);
    get_kb_i(p_im1_equil_norm, p_i_equil_norm, kb_i_equil);
    
    float weight = get_weights(p_im1, p_i);
    float equil_weight = get_weights(p_im1_equil, p_i_equil);
    
    float mutual_l[3];
    float equil_mutual_l[3];
    get_mutual_element_inverse(p_im1, p_i, weight, OUT mutual_l);
    get_mutual_element_inverse(p_im1_equil, p_i_equil, weight, OUT equil_mutual_l);
    
    float m_im1_transported[3];
    float m_im1_equil_transported[3];
    parallel_transport(m_im1, m_im1_transported, p_im1_norm, mutual_l);
    parallel_transport(m_im1_equil, m_im1_equil_transported, p_im1_equil_norm, equil_mutual_l);
    
    float m_i_transported[3];
    float m_i_equil_transported[3];
    parallel_transport(m_i, m_i_transported, p_i_norm, mutual_l);
    parallel_transport(m_i_equil, m_i_equil_transported, p_i_equil_norm, equil_mutual_l);
        
    float m_mutual[3];
    get_mutual_axes_inverse(m_im1_transported, m_i_transported, weight, m_mutual);

    float m_mutual_equil[3];
    get_mutual_axes_inverse(m_im1_equil_transported, m_i_equil_transported, equil_weight, m_mutual_equil);

    normalize(m_mutual_equil, m_mutual_equil);
    normalize(m_mutual, m_mutual);

    float n_mutual[3];
    float n_mutual_equil[3];    
    cross_product(mutual_l, m_mutual, n_mutual);
    cross_product(equil_mutual_l, m_mutual_equil, n_mutual_equil);
    
    // finally get omega
    //float omega_j_im1[2];    
    get_omega_j_i(kb_i, n_mutual, m_mutual, omega);
    
    //float omega_j_im1_equil[2];    
    get_omega_j_i(kb_i_equil, n_mutual_equil, m_mutual_equil, omega_equil);
    
    float bend_energy = 0;
    bend_energy += get_bend_energy(omega, omega_equil, B_i_equil);
    bend_energy = bend_energy*(0.5/(L_i)); // constant!

    return bend_energy;
}

int main(){
 return 0;  
}
