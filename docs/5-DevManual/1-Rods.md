FFEA_rod developer documentation {#rodsdev}
=================================================

A C++ implementation of FFEA_rod is provided within the FFEA software package. The relevant files are:

* rod_math_v9.cpp and rod_math_v9.h - the core math functions needed to compute energies and dynamics for the rod
* rod_structure.cpp and rod_structure.h - the data structures, I/O and various functions to do with setting rod positions, initialising rods and performing timesteps.
* dimensions.h - the FFEA header file specifying FFEA internal units.

In addition, FFEA_rod depends upon RngStream.h by Pierre L'Ecuyer and Richard Simard, OpenMP, and Boost. If you wish to build FFEA_rod outside of FFEA, the files listed above and these three libraries should be enough for FFEA_rod to compile.

## Creating rod objects

All rod functions and objects live inside the rod:: namespace. The FFEA_rod object and data structure stores information specific to a particular rod, but not global simulation parameters (e.g. an instance of RNGstreams or the temperature). To create a rod object, use

```cpp
rod::Rod* current_rod = new rod::Rod(filename, rod_id);
current_rod->load_header(filename);
current_rod->load_contents(filename);
current_rod->set_units();
```
For a rod to function properly, the following global simulation parameters must be set:
```cpp
current_rod->viscosity = stokes_visc; // solvent viscosity
current_rod->timestep = dt; // timestep
current_rod->kT = kT; // thermal energy (in FFEA units, see dimensions.h)
```
In FFEA, these quantities come from the global world parameters specified in the .ffea script file. In FFEA, multiple rods can be specified to run in a single simulation. They are parallel per node, not per rod.

## Running rod dynamics

To run a step of rod dynamics, use
```cpp
current_rod->do_timestep(rng);
```
Where 'rng' is an array of pointers to rngstream objects. This array should be at least as large as the environment variable 'OMP_NUM_THREADS'.

Finally, to write the current state of the rod to the rod trajectory file, use
```cpp
current_rod->write_frame_to_file();
```
By default, this will write the trajectory out to the same file that was read in. There is no differece between rod input and output files (other than the fact that the rod output file contains multiple frames). It is often convenient to make the rod write the trajectroy to another file, and keep our input file as-is:
```cpp
current_rod->change_filename("new_filename.rodtraj");
```
## Other rod member functions
Translate the rod:
```cpp
float translate[3] = {1.0,1.0,1.0};
current_rod->translate_rod(current_rod->current_r, translate);
current_rod->translate_rod(current_rod-equil_r, translate);
```

Rotate the rod:
```cpp
float euler_angles[3] = {1.0*M_PI,0.0,0.5*M_PI};
current_rod->rotate_rod(euler_angles);
```
## Miscellaneous info
* Unlike the rest of FFEA, rods always run with the same precision - single-precision floating point. You can safely cast these numbers to FFEA scalars, and vice-versa.
* The rods make use of 1-d arrays, 3 elements in length, to represent vectors, unlike the rest of FFEA, which uses a special 'vector' class. This is done for speed and memory efficiency. The internal rod state (e.g. current_rod->current_r) is also stored as a 1-D array, of length n*3, where n is the number of nodes. The 5-node 'window' needed to compute the energies ends up in a 2-d array, where the first index is node number.
* x, y, and z are named constants corresponding to the integers 0, 1, and 2. So, for example, declaring the variable float node[3] = {1,2,3}; would mean that node[y] would equal 2.
* This is taken even further by the vec3d macro, which simply loops over the x, y and z elements in an 3-element array containing 3-d  positional data.
* In 5-node slice needed to compute energies, as used in get_energies, the four elements in the 5-node are referred to by their indices im2, im1, i and ip1 (i-2, i-1, i and i+1). These indices are also named constants, so you can have (for example) p[im1][x].
* Do not make and delete lots of rod elements! They don't have a destructor! Well, they do, it's just that the destructor is  'when the program is closed, the memory is freed'.
* The math functions in mat_vec_functions are not reused.
* OpemMP acts upon the per-node loop of the energy calculations, but not the dynamics, resulting in a 4% serial fraction.
* Unlike rod_math and rod_structure, the contents of rod_blob_interface are deeply coupled to the FFEA algorithm and associated data structures. Please refer instead to the publication on this topic (which isn't out yet) or contact us.

## Rod-blob connections

Unfortunately, rod-blob connections are tightly coupled to the existing FFEA data structures and also have a lot of state, so be careful when creating or using them.

```cpp
rod::Rod_blob_interface* curr_rbi = new rod::Rod_blob_interface(connected_rod_ptr, connected_blob_ptr, ends_at_rod, to_index, from_index, blob_node_ids, rotation, node_weighting, order);
```
Rod-blob connections (called 'interfaces' in the code) are initialised after the rods and the blobs have been created. The pointers to the FFEA_rod and FFEA_blob objects are given as parameters. ends_at_rod is a bool, if true then the connection is blob_to_rod, if false it's rod_to_blob. to_index and from_index are integers, the indices of the rod nodes/blob elements in the connection, the blob node ids is an array of three integers (the nodes on the connected face), rotation is an array of floats specifying the euler angles, node_weighting is an array of floats specifying the position of the connection within the element, and order is initialisation order ID, another integer.

You can't use a rod-blob connection right after connecting it. First, in order of order, update the internal state of the connection. Updating the internal state copies the state of the tetrahedron from the blob to the rod. This copying is necessary because the state of the tetrahedron is actually different when calculating the energies in the rod-blob connection.

```cpp
curr_rbi->update_internal_state(true, true);
```
If, for any reason, the connection is repositioned, this function must run. That includes relative positioning and situating objects in the simulation box, which happens later in the init process.

Then, position the connection. If the connection goes blob-to-rod (e.g. if ends_at_rod==true), then
```cpp
curr_rbi->position_rod_from_blob(false);
curr_rbi->position_rod_from_blob(true);
```
Or, if ends_at_rod=false:
```cpp
curr_rbi->position_blob_from_rod();
```

Finally, set the initial values. This includes the initial positions of the connection node and element, and the equilibrium jacobian of the tetrahedron. This should be done right before the simulation starts running.
```cpp
curr_rbi->set_initial_values();
```
Luckily, actually running the dynamics of the connection is a lot easier. It doesn't matter to much *when* you do them, but by default, they happen just after the rod dynamics, the very last thing that's done in a timestep.
```cpp
curr_rbi->do_connection_timestep();
```

