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
 *      Skeleton.h
 */

#ifndef SKEL_H_INCLUDED
#define SKEL_H_INCLUDED

#include <stdio.h>
#include <vector>

#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "FFEA_return_codes.h"

/*
 * The "Joint" class
 */
class Joint {

	public:

		Joint();
		~Joint();

		/** The joint location, which is defined by an element centroid */
		tetra_element_linear *element;
		vector3 *pos;
};

/*
 * The "Bone" class
 */
class Bone {

	public:

		Bone();
		~Bone();

		/** Each bone points to 2 joints */
		Joint *joint[2];

		/** Each bone points to a list of nodes */
		vector<mesh_node*> linkedNodes;
};

/*
 * The "Skeleton" class
 */
class Skeleton {

	public:

	    /**
	     * Skeleton constructors:
	     */
	    Skeleton();
	    Skeleton(int num_joints, int num_bones);

	    /**
	     * Skeleton destructor:
	     */
	    ~Skeleton();

	    /**
	     * Method for adding a joint 
	     */
	    void add_joint(int j_index, int el_index, tetra_element_linear *elem);

	    /**
	     * Method for adding a bone 
	     */
	    void add_bone(int b_index, int j_index1, int j_index2);

	    /*
	     * Variables
	     */

	    /** Total number of bones (Skeleton elements) */
	    int num_bones;

	    /** Total number of joints (Skeleton nodes) */
	    int num_joints;

	    /** The bones, which simply point to the joints */
	    Bone *bone;

	    /** The joints, which are definined by element centroids */
	    Joint *joint;
};
#endif
