#ifndef VDW_SOLVER_H_INCLUDED
#define VDW_SOLVER_H_INCLUDED

#include <math.h>

#include "FFEA_return_codes.h"
#include "NearestNeighbourLinkedListCube.h"
#include "LJ_matrix.h"

class VdW_solver {

	public:
		VdW_solver()
		{
			total_num_surface_faces = 0;
			surface_face_lookup = NULL;
			box_size.x = 0;
			box_size.y = 0;
			box_size.z = 0;
		}
		~VdW_solver()
		{
			total_num_surface_faces = 0;
			surface_face_lookup = NULL;
			box_size.x = 0;
			box_size.y = 0;
			box_size.z = 0;
		}

		int init(NearestNeighbourLinkedListCube *surface_face_lookup, vector3 *box_size,  LJ_matrix *lj_matrix)
		{
			this->surface_face_lookup = surface_face_lookup;
			this->box_size.x = box_size->x;
			this->box_size.y = box_size->y;
			this->box_size.z = box_size->z;

			this->lj_matrix = lj_matrix;
			return FFEA_OK;
		}

		int solve()
		{
			const struct adjacent_cell_lookup_table_entry adjacent_cell_lookup_table[27] =
			{
				{-1, -1, -1},
				{-1, -1, 0},
				{-1, -1, +1},
				{-1, 0, -1},
				{-1, 0, 0},
				{-1, 0, +1},
				{-1, +1, -1},
				{-1, +1, 0},
				{-1, +1, +1},
				{0, -1, -1},
				{0, -1, 0},
				{0, -1, +1},
				{0, 0, -1},
				{0, 0, 0},
				{0, 0, +1},
				{0, +1, -1},
				{0, +1, 0},
				{0, +1, +1},
				{+1, -1, -1},
				{+1, -1, 0},
				{+1, -1, +1},
				{+1, 0, -1},
				{+1, 0, 0},
				{+1, 0, +1},
				{+1, +1, -1},
				{+1, +1, 0},
				{+1, +1, +1}
			};

			LinkedListNode<Face> *l_i = NULL;
			Face *f_i, *f_j;

			total_num_surface_faces = surface_face_lookup->get_pool_size();

			/* For each face, calculate the interaction with all other relevant faces and add the contribution to mat_C and mat_D */
			#ifdef FFEA_PER_BLOB_PARALLELISATION
				#pragma omp parallel for schedule(static)
			#endif
			for(int i = 0; i < total_num_surface_faces; i++) {

				// get the ith face
				l_i = surface_face_lookup->get_from_pool(i);
				f_i = l_i->obj;

				// Calculate this face's interaction with all faces in its cell and the 26 adjacent cells (3^3 = 27 cells)
				// Remember to check that the face is not interacting with itself or connected faces
				for(int c = 0; c < 27; c++) {
					LinkedListNode<Face> *l_j = NULL;
					l_j = surface_face_lookup->get_top_of_stack(	l_i->x + adjacent_cell_lookup_table[c].ix,
											l_i->y + adjacent_cell_lookup_table[c].iy,
											l_i->z + adjacent_cell_lookup_table[c].iz);
					while(l_j != NULL) {
						if(l_i->index != l_j->index) {
							f_j = l_j->obj;
							if(f_i->daddy_blob != f_j->daddy_blob) {
								//printf("(%d %d)\n", l_i->index, l_j->index);
								//if((l_i->index == 1 && l_j->index == 5) || (l_i->index == 5 && l_j->index == 1))
								do_interaction(f_i, f_j);
							}
						}
						l_j = l_j->next;
					}
				}

			}

			return FFEA_OK;
		}

		/* Allow protein VdW interactions along the top and bottom x-z planes */
		int solve_sticky_wall(scalar h)
		{
			int Nx = 0, Ny = 0, Nz = 0;
			surface_face_lookup->get_dim(&Nx, &Ny, &Nz);
			LinkedListNode<Face> *l_j = NULL;
			Face *f_j = NULL;
			for(int y = 0; y < Ny; y += Ny-1) {
				for(int z = 0; z < Nz; z++) {
					for(int x = 0; x < Nx; x++) {
						l_j = surface_face_lookup->get_top_of_stack(x, y, z);
						while(l_j != NULL) {
							f_j = l_j->obj;
							//fprintf(stderr, "Hi\n");
							do_sticky_xz_interaction(f_j, (y == 0), h * Ny);
							l_j = l_j->next;
						}
					}
				}
			}
			return FFEA_OK;
		}

		void get_interblob_vdw_stuff(Blob *dad_blob_1, Blob *dad_blob_2, scalar *total_area, scalar *total_force_mag, scalar *total_energy)
		{
			const struct adjacent_cell_lookup_table_entry adjacent_cell_lookup_table[27] =
			{
				{-1, -1, -1},
				{-1, -1, 0},
				{-1, -1, +1},
				{-1, 0, -1},
				{-1, 0, 0},
				{-1, 0, +1},
				{-1, +1, -1},
				{-1, +1, 0},
				{-1, +1, +1},
				{0, -1, -1},
				{0, -1, 0},
				{0, -1, +1},
				{0, 0, -1},
				{0, 0, 0},
				{0, 0, +1},
				{0, +1, -1},
				{0, +1, 0},
				{0, +1, +1},
				{+1, -1, -1},
				{+1, -1, 0},
				{+1, -1, +1},
				{+1, 0, -1},
				{+1, 0, 0},
				{+1, 0, +1},
				{+1, +1, -1},
				{+1, +1, 0},
				{+1, +1, +1}
			};
			LinkedListNode<Face> *l_i = NULL;
			Face *f_i, *f_j;

			total_num_surface_faces = surface_face_lookup->get_pool_size();
			vector3 total_force = {0.0, 0.0, 0.0};

			/* For each face, calculate the interaction with all other relevant faces and add the contribution to mat_C and mat_D */
			#ifdef FFEA_PER_BLOB_PARALLELISATION
				#pragma omp parallel for schedule(static)
			#endif
			for(int i = 0; i < total_num_surface_faces; i++) {

				// get the ith face
				l_i = surface_face_lookup->get_from_pool(i);
				f_i = l_i->obj;

				// Calculate this face's interaction with all faces in its cell and the 26 adjacent cells (3^3 = 27 cells)
				// Remember to check that the face is not interacting with itself or connected faces
				for(int c = 0; c < 27; c++) {
					LinkedListNode<Face> *l_j = NULL;
					l_j = surface_face_lookup->get_top_of_stack(	l_i->x + adjacent_cell_lookup_table[c].ix,
											l_i->y + adjacent_cell_lookup_table[c].iy,
											l_i->z + adjacent_cell_lookup_table[c].iz);
					while(l_j != NULL) {
						if(l_i->index != l_j->index) {
							f_j = l_j->obj;
							if((f_i->daddy_blob == dad_blob_1 && f_j->daddy_blob == dad_blob_2) || (f_i->daddy_blob == dad_blob_2 && f_j->daddy_blob == dad_blob_1)) {
					//			printf("(%d %d)\n", l_i->index, l_j->index);
								//if((l_i->index == 1 && l_j->index == 5) || (l_i->index == 5 && l_j->index == 1))
								do_measurement_interaction(f_i, f_j, total_area, &total_force, total_energy);
							}
						}
						l_j = l_j->next;
					}
				}

			}
			*total_force_mag = sqrt(total_force.x * total_force.x + total_force.y * total_force.y + total_force.z * total_force.z);

		}

	private:
		int total_num_surface_faces;
		NearestNeighbourLinkedListCube *surface_face_lookup;

		vector3 box_size;

		LJ_matrix *lj_matrix;

		struct adjacent_cell_lookup_table_entry {
			int ix, iy, iz;
		};

		struct tri_gauss_point
		{
			scalar W;
			scalar eta[3];
		};

		void do_interaction(Face *f1, Face *f2)
		{
			// First check two things (either of which results in not having to calculate anything):
			// Check that faces are facing each other, if not then they are not interacting
			if(dot(&f1->normal, &f2->normal) > 0) {
//				printf("DENIED NORMAL\n");
				return;
			}

			// Check that faces are in front of each other
//			vector3 sep = {f2->centroid.x - f1->centroid.x, f2->centroid.y - f1->centroid.y, f2->centroid.z - f1->centroid.z};
//			if(dot(&sep, &f1->normal) < 0 && dot(&sep, &f2->normal) > 0) {
//				printf("DENIED BEHIND\n");
//				return;
//			}

			// If faces are more than a few nanometres apart, don't bother calculating the force (it will be tiny)
//			if(distance2(&f1->centroid, &f2->centroid) > (vdw_r_eq * 6 * vdw_r_eq * 6)) {
//				return;
//			}
			
			// Get the interaction LJ parameters for these two face types
			scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
			lj_matrix->get_LJ_params(f1->vdw_interaction_type, f2->vdw_interaction_type, &vdw_eps, &vdw_r_eq);

			const int num_tri_gauss_quad_points = 3;

			const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] =
				{
					// Weight, eta1, eta2, eta3
				        {0.333333333333333,     {0.666666666666667, 0.166666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.666666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.166666666666667, 0.666666666666667}}

/*
				        {0.109951743655322,     {0.816847572980459, 0.091576213509771, 0.091576213509771}},
				        {0.109951743655322,     {0.091576213509771, 0.816847572980459, 0.091576213509771}},
				        {0.109951743655322,     {0.091576213509771, 0.091576213509771, 0.816847572980459}},
				        {0.223381589678011,     {0.108103018168070, 0.445948490915965, 0.445948490915965}},
				        {0.223381589678011,     {0.445948490915965, 0.108103018168070, 0.445948490915965}},
				        {0.223381589678011,     {0.445948490915965, 0.445948490915965, 0.108103018168070}}
*/

/*
					{0.050844906370207,     {0.873821971016996, 0.063089014491502, 0.063089014491502}},
					{0.050844906370207,     {0.063089014491502, 0.873821971016996, 0.063089014491502}},
					{0.050844906370207,     {0.063089014491502, 0.063089014491502, 0.873821971016996}},
					{0.116786275726379,     {0.501426509658179, 0.249286745170910, 0.249286745170910}},
					{0.116786275726379,     {0.249286745170910, 0.501426509658179, 0.249286745170910}},
					{0.116786275726379,     {0.249286745170910, 0.249286745170910, 0.501426509658179}},
					{0.082851075618374,     {0.636502499121399, 0.310352451033785, 0.053145049844816}},
					{0.082851075618374,     {0.310352451033785, 0.053145049844816, 0.636502499121399}},
					{0.082851075618374,     {0.053145049844816, 0.636502499121399, 0.310352451033785}},
					{0.082851075618374,     {0.636502499121399, 0.053145049844816, 0.310352451033785}},
					{0.082851075618374,     {0.310352451033785, 0.636502499121399, 0.053145049844816}},
					{0.082851075618374,     {0.053145049844816, 0.310352451033785, 0.636502499121399}}
*/
				};

			vector3 p[num_tri_gauss_quad_points], q[num_tri_gauss_quad_points];
			vector3 force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

			// Convert all area coordinate gauss points to cartesian
			for(int i = 0; i < num_tri_gauss_quad_points; i++) {
				f1->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
				f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);
			}

			// Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
				for(int l = k; l < num_tri_gauss_quad_points; l++) {
//					vector3 r = 	{
//								minimum_image(p[k].x - q[l].x, box_size.x),
//								minimum_image(p[k].y - q[l].y, box_size.y),
//								minimum_image(p[k].z - q[l].z, box_size.z)
//							};
					vector3 r = 	{
								p[k].x - q[l].x,
								p[k].y - q[l].y,
								p[k].z - q[l].z
							};


					scalar mag_r = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					scalar force_mag = 12 * pow(vdw_r_eq, 6) * vdw_eps * (pow(mag_r, -7) - pow(vdw_r_eq, 6) * pow(mag_r, -13));
					force_mag *= -1;

					force_pair_matrix[k][l].x = force_mag * (r.x / mag_r);
					force_pair_matrix[k][l].y = force_mag * (r.y / mag_r);
					force_pair_matrix[k][l].z = force_mag * (r.z / mag_r);

					force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
					force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
					force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
				}
			}

//			printf("YO\n");
//			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
//				for(int l = 0; l < num_tri_gauss_quad_points; l++) {
//					printf("(%e %e %e) ", force_pair_matrix[k][l].x, force_pair_matrix[k][l].y, force_pair_matrix[k][l].z);
//				}
//				printf("\n");
//			}

			scalar ApAq = f1->area * f2->area;

			for(int j = 0; j < 3; j++) {
				vector3 force1 = {0, 0, 0}, force2 = {0, 0, 0};
				for(int k = 0; k < num_tri_gauss_quad_points; k++) {
					for(int l = 0; l < num_tri_gauss_quad_points; l++) {
						scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
						scalar d = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];
						//printf("c = %e, %e, %e, %e\n", c, gauss_points[k].W, gauss_points[l].W, gauss_points[l].eta[j]);
						force1.x += c * force_pair_matrix[k][l].x;
						force1.y += c * force_pair_matrix[k][l].y;
						force1.z += c * force_pair_matrix[k][l].z;

						force2.x -= d * force_pair_matrix[l][k].x;
						force2.y -= d * force_pair_matrix[l][k].y;
						force2.z -= d * force_pair_matrix[l][k].z;
					}
				}
				force1.x *= ApAq;
				force1.y *= ApAq;
				force1.z *= ApAq;
				f1->add_force_to_node(j, &force1);

//				printf("1:: %d %e %e %e\n", j, force1.x, force1.y, force1.z);

				force2.x *= ApAq;
				force2.y *= ApAq;
				force2.z *= ApAq;
				f2->add_force_to_node(j, &force2);
//				printf("2:: %d %e %e %e\n", j, force2.x, force2.y, force2.z);

//				f2->add_force_to_node_atomic(j, &force);

//				if(j == 0) {
//					f2->add_force_to_node_atomic(0, &force);
//				} else if(j == 1) {
//					f2->add_force_to_node_atomic(1, &force);
//				} else if(j == 2) {
//					f2->add_force_to_node_atomic(2, &force);
//				} else {
//					printf("WTF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
//				}

//				printf("Face 1: force on node %d = {%e %e %e}\n", j, force.x, force.y, force.z);
			}
			

//			printf("\n\n\n\n");

/*
			for(int j = 0; j < 3; j++) {
				vector3 force = {0, 0, 0};
				for(int k = 0; k < num_tri_gauss_quad_points; k++) {
					for(int l = 0; l < num_tri_gauss_quad_points; l++) {
						scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];
						force.x -= c * force_pair_matrix[k][l].x;
						force.y -= c * force_pair_matrix[k][l].y;
						force.z -= c * force_pair_matrix[k][l].z;
					}
				}
				force.x *= ApAq;
				force.y *= ApAq;
				force.z *= ApAq;
				f2->add_force_to_node(j, &force);
				printf("Face 2: force on node %d = {%e %e %e}\n", j, force.x, force.y, force.z);
			}
			printf("\n");

*/					

		}

		void do_measurement_interaction(Face *f1, Face *f2, scalar *total_area, vector3 *total_force, scalar *total_energy)
		{
			// First check two things (either of which results in not having to calculate anything):
			// Check that faces are facing each other, if not then they are not interacting
			if(dot(&f1->normal, &f2->normal) > 0) {
//				printf("DENIED NORMAL\n");
				return;
			}

			// Check that faces are in front of each other
//			vector3 sep = {f2->centroid.x - f1->centroid.x, f2->centroid.y - f1->centroid.y, f2->centroid.z - f1->centroid.z};
//			if(dot(&sep, &f1->normal) < 0 && dot(&sep, &f2->normal) > 0) {
//				printf("DENIED BEHIND\n");
//				return;
//			}

			// If faces are more than a few nanometres apart, don't bother calculating the force (it will be tiny)
//			if(distance2(&f1->centroid, &f2->centroid) > (vdw_r_eq * 6 * vdw_r_eq * 6)) {
//				return;
//			}

			// Get the interaction LJ parameters for these two face types
			scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
			lj_matrix->get_LJ_params(f1->vdw_interaction_type, f2->vdw_interaction_type, &vdw_eps, &vdw_r_eq);

			const int num_tri_gauss_quad_points = 3;

			const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] =
				{
					// Weight, eta1, eta2, eta3
				        {0.333333333333333,     {0.666666666666667, 0.166666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.666666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.166666666666667, 0.666666666666667}}

/*
				        {0.109951743655322,     {0.816847572980459, 0.091576213509771, 0.091576213509771}},
				        {0.109951743655322,     {0.091576213509771, 0.816847572980459, 0.091576213509771}},
				        {0.109951743655322,     {0.091576213509771, 0.091576213509771, 0.816847572980459}},
				        {0.223381589678011,     {0.108103018168070, 0.445948490915965, 0.445948490915965}},
				        {0.223381589678011,     {0.445948490915965, 0.108103018168070, 0.445948490915965}},
				        {0.223381589678011,     {0.445948490915965, 0.445948490915965, 0.108103018168070}}
*/

/*
					{0.050844906370207,     {0.873821971016996, 0.063089014491502, 0.063089014491502}},
					{0.050844906370207,     {0.063089014491502, 0.873821971016996, 0.063089014491502}},
					{0.050844906370207,     {0.063089014491502, 0.063089014491502, 0.873821971016996}},
					{0.116786275726379,     {0.501426509658179, 0.249286745170910, 0.249286745170910}},
					{0.116786275726379,     {0.249286745170910, 0.501426509658179, 0.249286745170910}},
					{0.116786275726379,     {0.249286745170910, 0.249286745170910, 0.501426509658179}},
					{0.082851075618374,     {0.636502499121399, 0.310352451033785, 0.053145049844816}},
					{0.082851075618374,     {0.310352451033785, 0.053145049844816, 0.636502499121399}},
					{0.082851075618374,     {0.053145049844816, 0.636502499121399, 0.310352451033785}},
					{0.082851075618374,     {0.636502499121399, 0.053145049844816, 0.310352451033785}},
					{0.082851075618374,     {0.310352451033785, 0.636502499121399, 0.053145049844816}},
					{0.082851075618374,     {0.053145049844816, 0.310352451033785, 0.636502499121399}}
*/
				};

			vector3 p[num_tri_gauss_quad_points], q[num_tri_gauss_quad_points];
			vector3 force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

			// Convert all area coordinate gauss points to cartesian
			for(int i = 0; i < num_tri_gauss_quad_points; i++) {
				f1->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
				f2->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &q[i]);
			}
		
			// Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
				for(int l = k; l < num_tri_gauss_quad_points; l++) {
//					vector3 r = 	{
//								minimum_image(p[k].x - q[l].x, box_size.x),
//								minimum_image(p[k].y - q[l].y, box_size.y),
//								minimum_image(p[k].z - q[l].z, box_size.z)
//							};
					vector3 r = 	{
								p[k].x - q[l].x,
								p[k].y - q[l].y,
								p[k].z - q[l].z
							};


					scalar mag_r = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
					scalar force_mag = 12 * pow(vdw_r_eq, 6) * vdw_eps * (pow(mag_r, -7) - pow(vdw_r_eq, 6) * pow(mag_r, -13));
					
					*total_energy += vdw_eps * pow(vdw_r_eq/mag_r, 6) * (pow(vdw_r_eq/mag_r, 6) - 2);
					force_mag *= -1;

					force_pair_matrix[k][l].x = force_mag * (r.x / mag_r);
					force_pair_matrix[k][l].y = force_mag * (r.y / mag_r);
					force_pair_matrix[k][l].z = force_mag * (r.z / mag_r);

					force_pair_matrix[l][k].x = force_pair_matrix[k][l].x;
					force_pair_matrix[l][k].y = force_pair_matrix[k][l].y;
					force_pair_matrix[l][k].z = force_pair_matrix[k][l].z;
				}
			}

			scalar ApAq = f1->area * f2->area;
			*total_area += f1->area + f2->area;
			*total_energy *= f1->area * f2->area;

			for(int j = 0; j < 3; j++) {
				vector3 force1 = {0, 0, 0}, force2 = {0, 0, 0};
				for(int k = 0; k < num_tri_gauss_quad_points; k++) {
					for(int l = 0; l < num_tri_gauss_quad_points; l++) {
						scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
						scalar d = gauss_points[k].W * gauss_points[l].W * gauss_points[k].eta[j];
						//printf("c = %e, %e, %e, %e\n", c, gauss_points[k].W, gauss_points[l].W, gauss_points[l].eta[j]);
						force1.x += c * force_pair_matrix[k][l].x;
						force1.y += c * force_pair_matrix[k][l].y;
						force1.z += c * force_pair_matrix[k][l].z;

						force2.x -= d * force_pair_matrix[l][k].x;
						force2.y -= d * force_pair_matrix[l][k].y;
						force2.z -= d * force_pair_matrix[l][k].z;
					}
				}
				force2.x *= ApAq;
				force2.y *= ApAq;
				force2.z *= ApAq;
				total_force->x += force2.x;
				total_force->y += force2.y;
				total_force->z += force2.z;
				
			} 
				

		}

		void do_sticky_xz_interaction(Face *f, bool bottom_wall, scalar dim_y)
		{
			scalar y_wall = 0;//-vdw_r_eq;
			if(bottom_wall == false) {
				y_wall = dim_y;// + vdw_r_eq;
			}

			// Check that face is facing wall, if not then it should not be interacting
			if((bottom_wall == true && f->normal.y > 0) || (bottom_wall == false && f->normal.y < 0)) {
				return;
			}

			// Get the interaction LJ parameters for these two face types
			scalar vdw_eps = 0.0, vdw_r_eq = 0.0;
			lj_matrix->get_LJ_params(f->vdw_interaction_type, f->vdw_interaction_type, &vdw_eps, &vdw_r_eq);

			const int num_tri_gauss_quad_points = 3;

			const struct tri_gauss_point gauss_points[num_tri_gauss_quad_points] =
				{
					// Weight, eta1, eta2, eta3
				        {0.333333333333333,     {0.666666666666667, 0.166666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.666666666666667, 0.166666666666667}},
				        {0.333333333333333,     {0.166666666666667, 0.166666666666667, 0.666666666666667}}

/*
				        {0.109951743655322,     {0.816847572980459, 0.091576213509771, 0.091576213509771}},
				        {0.109951743655322,     {0.091576213509771, 0.816847572980459, 0.091576213509771}},
				        {0.109951743655322,     {0.091576213509771, 0.091576213509771, 0.816847572980459}},
				        {0.223381589678011,     {0.108103018168070, 0.445948490915965, 0.445948490915965}},
				        {0.223381589678011,     {0.445948490915965, 0.108103018168070, 0.445948490915965}},
				        {0.223381589678011,     {0.445948490915965, 0.445948490915965, 0.108103018168070}}
*/

/*
					{0.050844906370207,     {0.873821971016996, 0.063089014491502, 0.063089014491502}},
					{0.050844906370207,     {0.063089014491502, 0.873821971016996, 0.063089014491502}},
					{0.050844906370207,     {0.063089014491502, 0.063089014491502, 0.873821971016996}},
					{0.116786275726379,     {0.501426509658179, 0.249286745170910, 0.249286745170910}},
					{0.116786275726379,     {0.249286745170910, 0.501426509658179, 0.249286745170910}},
					{0.116786275726379,     {0.249286745170910, 0.249286745170910, 0.501426509658179}},
					{0.082851075618374,     {0.636502499121399, 0.310352451033785, 0.053145049844816}},
					{0.082851075618374,     {0.310352451033785, 0.053145049844816, 0.636502499121399}},
					{0.082851075618374,     {0.053145049844816, 0.636502499121399, 0.310352451033785}},
					{0.082851075618374,     {0.636502499121399, 0.053145049844816, 0.310352451033785}},
					{0.082851075618374,     {0.310352451033785, 0.636502499121399, 0.053145049844816}},
					{0.082851075618374,     {0.053145049844816, 0.310352451033785, 0.636502499121399}}
*/
				};

			vector3 p[num_tri_gauss_quad_points];
			scalar force_pair_matrix[num_tri_gauss_quad_points][num_tri_gauss_quad_points];

			// Convert all area coordinate gauss points to cartesian
			for(int i = 0; i < num_tri_gauss_quad_points; i++) {
				f->barycentric_calc_point(gauss_points[i].eta[0], gauss_points[i].eta[1], gauss_points[i].eta[2], &p[i]);
			}

			// Construct the force pair matrix: f(p, q) where p and q are all the gauss points in each face
			for(int k = 0; k < num_tri_gauss_quad_points; k++) {
				for(int l = k; l < num_tri_gauss_quad_points; l++) {
					scalar mag_r = p[k].y - y_wall;

					scalar force_mag = 12 * pow(vdw_r_eq, 6) * vdw_eps * (pow(mag_r, -7) - pow(vdw_r_eq, 6) * pow(mag_r, -13));
					force_mag *= -1;

					force_pair_matrix[k][l] = force_mag;
					force_pair_matrix[l][k] = force_pair_matrix[k][l];
				}
			}

			scalar Asq = f->area * f->area;
			for(int j = 0; j < 3; j++) {
				vector3 force = {0, 0, 0};
				for(int k = 0; k < num_tri_gauss_quad_points; k++) {
					for(int l = 0; l < num_tri_gauss_quad_points; l++) {
						scalar c = gauss_points[k].W * gauss_points[l].W * gauss_points[l].eta[j];
						force.y += c * force_pair_matrix[k][l];
					}
				}
				force.y *= Asq;
				f->add_force_to_node(j, &force);
			}

			// set vdw interaction flag for this face
			f->set_vdw_xz_interaction_flag(true);
		}


		scalar distance2(vector3 *p, vector3 *q)
		{
			scalar dx = p->x - q->x, dy = p->y - q->y, dz = p->z - q->z;
			return dx * dx + dy * dy + dz * dz;
		}

		scalar dot(vector3 *p, vector3 *q)
		{
			return p->x * q->x + p->y * q->y + p->z * q->z;
		}

		scalar dot_with_normal(vector3 *p, vector3 *q, vector3 *n)
		{
			return (q->x - p->x) * n->x + (q->y - p->y) * n->y + (q->z - p->z) * n->z;
		}

		scalar minimum_image(scalar delta, scalar size)
		{
			if(fabs(delta) > size * .5) {
				return size - delta;
			}

			return delta;
		}
};

#endif
