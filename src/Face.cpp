#include "Face.h"

Face::Face() {
    n[0] = NULL;
    n[1] = NULL;
    n[2] = NULL;
    n[3] = NULL;
    e = NULL;
    vdw_interaction_type = -1;
    area_0 = 0;
    zero_force();
    num_blobs = 0;
    vdw_xz_interaction_flag = false;
    vdw_bb_interaction_flag = NULL;
    vdw_bb_force = NULL;
    vdw_bb_energy = NULL;
    vdw_xz_force = NULL;
    vdw_xz_energy = 0.0;
    daddy_blob = NULL;
}

Face::~Face() {
    n[0] = NULL;
    n[1] = NULL;
    n[2] = NULL;
    n[3] = NULL; 
    e = NULL;
    vdw_interaction_type = -1;
    area_0 = 0;
    zero_force();
    num_blobs = 0;
    vdw_xz_interaction_flag = false;
    delete[] vdw_bb_interaction_flag;
    vdw_bb_interaction_flag = NULL;

    delete[] vdw_bb_force;
    vdw_bb_force = NULL;
    delete[] vdw_bb_energy;
    vdw_bb_energy = NULL;
    delete[] vdw_xz_force;
    vdw_xz_force = NULL;
    vdw_xz_energy = 0.0;
    daddy_blob = NULL;
}

void Face::init(tetra_element_linear *e, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *opposite, SecondOrderFunctions::stu centroid_stu, Blob *daddy_blob, SimulationParams *params) {
    this->e = e;
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
    n[3] = opposite; 

    calc_area_normal_centroid();
    area_0 = area;

    this->centroid_stu.s = centroid_stu.s;
    this->centroid_stu.t = centroid_stu.t;
    this->centroid_stu.u = centroid_stu.u;

    this->num_blobs = params->num_blobs;
    vdw_bb_force = new vector3[num_blobs];
    vdw_bb_energy = new scalar[num_blobs];
    vdw_bb_interaction_flag = new bool[num_blobs];
    vdw_xz_force = new vector3;

    for(int i = 0; i < num_blobs; ++i) {
	vector3_set_zero(&vdw_bb_force[i]);
        vdw_bb_energy[i] = 0.0;
	vdw_bb_interaction_flag[i] = false;
    }
    vector3_set_zero(vdw_xz_force);
    vdw_xz_energy = 0.0;


    this->daddy_blob = daddy_blob;
}

void Face::init(mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *opposite, Blob *daddy_blob, SimulationParams *params) {
    this->e = NULL;
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
    n[3] = opposite;

    calc_area_normal_centroid();
    area_0 = area;

    this->centroid_stu.s = 0;
    this->centroid_stu.t = 0;
    this->centroid_stu.u = 0;

    this->num_blobs = params->num_blobs;
    vdw_bb_force = new vector3[num_blobs];
    vdw_bb_energy = new scalar[num_blobs];
    vdw_bb_interaction_flag = new bool[num_blobs];
    vdw_xz_force = new vector3;

    for(int i = 0; i < num_blobs; ++i) {
	vector3_set_zero(&vdw_bb_force[i]);
        vdw_bb_energy[i] = 0.0;
	vdw_bb_interaction_flag[i] = false;
    }
    vector3_set_zero(vdw_xz_force);
    vdw_xz_energy = 0.0;

    this->daddy_blob = daddy_blob;
}

void Face::set_vdw_interaction_type(int vdw_interaction_type) {
    this->vdw_interaction_type = vdw_interaction_type;
}

void Face::calc_area_normal_centroid() {
    // (1/2) * |a x b|
    vector3 a = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z},
    b = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};
    normal.x = a.y * b.z - a.z * b.y;
    normal.y = a.z * b.x - a.x * b.z;
    normal.z = a.x * b.y - a.y * b.x;

    scalar normal_mag = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

    area = .5 * normal_mag;

    // Normalise normal
    normal.x /= normal_mag;
    normal.y /= normal_mag;
    normal.z /= normal_mag;

    // Find centroid
    centroid.x = (1.0 / 3.0) * (n[0]->pos.x + n[1]->pos.x + n[2]->pos.x);
    centroid.y = (1.0 / 3.0) * (n[0]->pos.y + n[1]->pos.y + n[2]->pos.y);
    centroid.z = (1.0 / 3.0) * (n[0]->pos.z + n[1]->pos.z + n[2]->pos.z);
}

vector3 * Face::get_centroid() {
	centroid.x = (1.0 / 3.0) * (n[0]->pos.x + n[1]->pos.x + n[2]->pos.x);
	centroid.y = (1.0 / 3.0) * (n[0]->pos.y + n[1]->pos.y + n[2]->pos.y);
	centroid.z = (1.0 / 3.0) * (n[0]->pos.z + n[1]->pos.z + n[2]->pos.z);
	
	return &centroid;
}

scalar Face::get_area() {

    vector3 temp;

    // (1/2) * |a x b|
    vector3 a = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z},
    b = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};
    temp.x = a.y * b.z - a.z * b.y;
    temp.y = a.z * b.x - a.x * b.z;
    temp.z = a.x * b.y - a.y * b.x;

    scalar normal_mag = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

    area = .5 * normal_mag;
    return area;

}

/* Calculate the point p on this triangle given the barycentric coordinates b1, b2, b3 */
void Face::barycentric_calc_point(scalar b1, scalar b2, scalar b3, vector3 *p) {
    p->x = b1 * n[0]->pos.x + b2 * n[1]->pos.x + b3 * n[2]->pos.x;
    p->y = b1 * n[0]->pos.y + b2 * n[1]->pos.y + b3 * n[2]->pos.y;
    p->z = b1 * n[0]->pos.z + b2 * n[1]->pos.z + b3 * n[2]->pos.z;
}

/* Returns the average electrostatic potential of this face */
scalar Face::average_phi() {
    return (1.0 / 3.0) * (n[0]->phi + n[1]->phi + n[2]->phi);
}

scalar Face::get_normal_flux() {
    vector3 dphi;
    e->get_grad_phi_at_stu(&dphi, centroid_stu.s, centroid_stu.t, centroid_stu.u);
    return dphi.x * normal.x + dphi.y * normal.y + dphi.z * normal.z;
}

/*
                scalar get_normal_flux()
                {
                        vector3 dphi =  {
                                        e->n[0]->phi * e->dpsi[0] + e->n[1]->phi * e->dpsi[1] + e->n[2]->phi * e->dpsi[2] + e->n[3]->phi * e->dpsi[3],
                                        e->n[0]->phi * e->dpsi[4] + e->n[1]->phi * e->dpsi[5] + e->n[2]->phi * e->dpsi[6] + e->n[3]->phi * e->dpsi[7],
                                        e->n[0]->phi * e->dpsi[8] + e->n[1]->phi * e->dpsi[9] + e->n[2]->phi * e->dpsi[10] + e->n[3]->phi * e->dpsi[11]
                                        };
                        return dphi.x * normal.x + dphi.y * normal.y + dphi.z * normal.z;
                }
 */


void Face::add_force_to_node(int i, vector3 *f) {
    force[i].x += f->x;
    force[i].y += f->y;
    force[i].z += f->z;
}

void Face::add_force_to_node_atomic(int i, vector3 *f) {
    force[i].x += f->x;
    force[i].y += f->y;
    force[i].z += f->z;
}

void Face::add_bb_vdw_force_to_record(vector3 *f, int other_blob_index) {
    vdw_bb_force[other_blob_index].x += f->x;
    vdw_bb_force[other_blob_index].y += f->y;
    vdw_bb_force[other_blob_index].z += f->z;
}

void Face::add_bb_vdw_energy_to_record(scalar energy, int other_blob_index) {
    vdw_bb_energy[other_blob_index] += energy;
}

void Face::add_xz_vdw_force_to_record(vector3 *f) {
    vdw_xz_force->x += f->x;
    vdw_xz_force->y += f->y;
    vdw_xz_force->z += f->z;
}

void Face::add_xz_vdw_energy_to_record(scalar energy) {
    vdw_xz_energy += energy;
}

void Face::zero_force() {
    for (int i = 0; i < 3; i++) {
        force[i].x = 0;
        force[i].y = 0;
        force[i].z = 0;
    }

}

void Face::zero_vdw_bb_measurement_data() {
    for (int i = 0; i < num_blobs; ++i) {
        vdw_bb_force[i].x = 0.0;
        vdw_bb_force[i].y = 0.0;
        vdw_bb_force[i].z = 0.0;
        vdw_bb_energy[i] = 0.0;
    }
}

void Face::zero_vdw_xz_measurement_data() {
    vdw_xz_force->x = 0.0;
    vdw_xz_force->y = 0.0;
    vdw_xz_force->z = 0.0;
    vdw_xz_energy = 0.0;
}

void Face::set_vdw_xz_interaction_flag(bool state) {
    vdw_xz_interaction_flag = state;
}

void Face::set_vdw_bb_interaction_flag(bool state, int other_blob_index) {
    vdw_bb_interaction_flag[other_blob_index] = state;
}

bool Face::is_vdw_active() {
    if (vdw_interaction_type == -1) {
        return false;
    } else {
        return true;
    }
}

bool Face::checkTetraIntersection(Face *f2) {
  
  // V1 = n
  // V2 = f2->n
  scalar tetA[4][3], tetB[4][3]; 
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  return (tet_a_tet(tetA, tetB)); 
 
} 

scalar Face::getTetraIntersectionVolume(Face *f2){

  scalar tetA[4][3], tetB[4][3]; 
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  return volumeIntersection(tetA, tetB);


} 
