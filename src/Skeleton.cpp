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

#include "Skeleton.h"

Skeleton::Skeleton() {
	num_joints = 0;
	num_bones = 0;
	joint = NULL;
	bone = NULL;
 
}

Skeleton::Skeleton(int num_joints, int num_bones) {
	this->num_joints = num_joints;
	this->num_bones = num_bones;
	joint = new Joint[num_joints];
	bone = new Bone[num_joints];
 
}

Skeleton::~Skeleton() {
	num_joints = 0;
	num_bones = 0;

	delete[] joint;
	joint = NULL;

	delete[] bone;
	bone = NULL;
}

void Skeleton::add_joint(int j_index, int el_index, tetra_element_linear *elem) {

	joint[j_index].element = &elem[el_index];
	joint[j_index].pos = &joint[j_index].element->centroid;
}

void Skeleton::add_bone(int b_index, int j_index1, int j_index2) {

	bone[b_index].joint[0] = &joint[j_index1];
	bone[b_index].joint[1] = &joint[j_index2];
	bone[b_index].calculate_centroid();
}

Joint::Joint() {
	element = NULL;
	pos = NULL;
}

Joint::~Joint() {
	element = NULL;
	pos = NULL;
}

Bone::Bone() {

	for(int i = 0; i < 2; ++i) {
		joint[i] = NULL;
	}
	linkedNodes.clear();
	vector3_set_zero(centroid);
}

Bone::~Bone() {

	for(int i = 0; i < 2; ++i) {
		joint[i] = NULL;
	}
	linkedNodes.clear();
	vector3_set_zero(centroid);
}

void Bone::calculate_centroid() {
	for(int i = 0; i < 3; ++i) {
		centroid[i] = 0.5 * ((*joint[0]->pos)[i] + (*joint[1]->pos)[i]);
	}
}
