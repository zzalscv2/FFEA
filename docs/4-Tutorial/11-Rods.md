Using FFEA_rod {#rods}
=========================

## The FFEA_rod object

You can create rods using the FFEA_rod python API. With FFEA installed, you can import it with

```python
import FFEA_rod # FFEA installed normally
import ffeatools.modules.FFEA_rod # FFEA python paackage installed
```
In order to manually create a rod, you can initialize a rod object without a filename:
```python
my_rod = FFEA_rod.FFEA_rod(num_elements=10) 
```
The attributes of the rod are stored in numpy `ndarray` objects. There are 3 dimensions, the frame, the element\node index, and the dimension. Rods created like this only have one frame. The relevant attributes are
*  current_r - current position of the nodes in the rod.
* current_m - current vectors for the material frames.
* equil_r - equilibrium configuration of the rod.
* equil_m - equilibrium material frame of the rod.
* perturbed x, y, z and twist energy, positive and negative - results from each stage of the numerical derivatives taken during the simulation. These are used for analysis.
* material_params - each node has the structure [stretch, twist, radius]. This mixes indices, as the stretch and radii are properties of elements, whereas the twist is the property of a node.
* B_matrix - the parameter for the bending energy. Each element of each frame is a 4-element array specifying the contents of the 2x2 B matrix (the bending modulus).

Try accessing and modifying some of these now! To write out a newly-created rod,
```python
my_rod.write_rod("my_rod.rod")
```
The FFEA_rod data format is actually a csv file which can be opened in a text editor. Try doing this now! By default, the rod simulation writes out files with the .rodtraj format - these are identical in format, they just contain more than one frame. As rod files are actually rod trajectories containing only one frame, when modifying the contents of the rod, we need to add a zero index to the rod attributes (e.g. `material_params[0]`).

# Creating rods

Although you can create a rod by manually by populating the numpy arrays listed above, the `FFEA_rod.rod_creator` object provides an easy way to populate some of these.

To create a rod from a PDB file, use `parse_pdb`.
```python
my_nodes = FFEA_rod.rod_creator.parse_pdb("my_pdb.pdb")
```
To select only particular atoms, provide a list of their names in the parameter `atom_names`. By default, the value of names is `['CA']` so that the backbone is selected. To take a look at the rod, use `FFEA_rod.rod_creator.preview_rod(my_nodes)`.

This will not create a rod on its own! This method will create one rod node for every atom in the PDB file that matches the `atom_names` filter. The data is output as a 2-d numpy array, representing a single frame of a rod trajectory. To coarse-grain, use create_rod_spline. This will automatically interpolate the points in the protein backbone as a b-spline.

```python
my_nodes_interpolated = FFEA_rod.rod_creator.create_rod_spline(my_nodes, 10)
```
This will create another 2D numpy array. If you know about b-splines, you can use the parameters `degree` and `smoothness`to further control the spline that comes out. You can again preview using `FFEA_rod.rod_creator.preview_rod(my_nodes_interpolated)`.

To set the contents of the rod, set the atrributes of the blank rod object you created earlier.

```python
my_rod.current_r[0] = my_nodes_interpolated
my_rod.equil_r[0] = my_nodes_interpolated # assuming the PDB is of the equiulibrium structure!
```

You can create the same kind of arrays from parametric functions. For example:
```python
import numpy

def x_func(t):
    return np.sin(t)
    
def y_func(t):
    return np.cos(t)
    
def z_func(t):
    return t
    
my_nodes = create_rod_parametric(x_func, y_func, z_func, 0, 1e-8, 10)
```
This creates a spiralling rod with 10 nodes which is 10nm in length.

Two functions are provided to set up material axes. To create the initial set of material axes, use `create_material_frame`. Note that this requires a rod object as a parameter, not an array of nodes.

```
my_rod.current_m[0], my_rod.equil_m[0] = FFEA_rod.rod_creator.create_material_frame(my_rod)
```
This will create a set of arbitrary material axis with no equilibrium. If, for some crazy reason, you want an equilibrium twist, you can add it like this:

```python
def rotate_function(t):
    return sin(t)
    
my_rod.current_m[0] = FFEA_rod.create_rod.rotate_material_frame(my_rod, rotate_function, my_rod.current_m, function_range(0, 2*pi))
```

To set material parameters on the rod, use `set_params`. For example, with a known stretch constant, torsion constant, radius and isotropic B matrix:
```python
# some sensible defaults
stretch_constant = 3.5e-11
twist_constant = 5e-29
radius = 5e-9
bend_constant = 3.5e-29
set_params(my_rod, stretch_constant, torsion_constant, radius, bending_modulus=bend_constant)
```
`set_params` also accepts the values of the young's modulus and radius to set the bending modulus, and you can specify a list of indices to apply these values to with the parameter `rod_segments`. Note again that the indexing will be off-by-one for properties specified at elements rather than nodes (such as the stretching constant). Sorry.

To write your finished rod to a file, run
```python
my_rod.write_rod("my_rod.rod")
```

## Adding a rod to an FFEA script
Rods can be included in FFEA script files either with or without other FFEA objects (although the two cannot yet interact). Here's an example of an FFEA script containing a rod:
```xml
<param>
	<restart = 0>
	<dt = 1e-11>
	<kT = 4.11e-21>
	<check = 1000>
	<num_steps = 20000000> 
	<rng_seed = time>
	<trajectory_out_fname = delet.ftj>
	<measurement_out_fname = delet.fm>
	<epsilon = 0.01>
	<max_iterations_cg = 1000>
	<kappa = 2e9>
	<epsilon_0 = 1>
	<dielec_ext = 1>
	<calc_stokes = 1>
	<stokes_visc = 1e-03>
	<calc_vdw = 1>
	<calc_noise = 1>
	<calc_kinetics = 0>
	<vdw_type = steric>
	<vdw_cutoff = 5e-10>
   <inc_self_vdw = 0>
	<vdw_steric_factor = 2e3>
	<calc_es = 0>
	<es_update = 4>
	<es_N_x = 50>
	<es_N_y = 30>
	<es_N_z = 30>
	<sticky_wall_xz = 0>
	<wall_x_1 = PBC>
	<wall_x_2 = PBC>
	<wall_y_1 = PBC>
	<wall_y_2 = PBC>
	<wall_z_1 = PBC>
	<wall_z_2 = PBC>
	<es_h = 1>
	<num_blobs = 0>
	<num_rods = 1>
	<num_conformations = (1)>
	<num_states = (1)>
</param>
<system>

    <rod>
    		<input = dest.rod>
    		<output = dest.rodtraj>
    		<centroid_pos = (0.0,0.0,0.0)> <!-- Translation vector -->
    		<rotation = (0.00, 0.00, 0,00)> <!-- Euler angles, X, Y and Z -->
    		<scale = 1>
    </rod>

</system>
```
Note: some of these parameters (such as restart) do not yet have an affect on rods. Additionally, rod simulations running inside FFEA will experience reduced performance compared to rod simulations running on their own. For information on how to set up an isolated rod simulation, please check the developer documentation.

Run a rod simulation the exact same way that you would run a regular FFEA simulation: e.g.
```bash
ffea ffeascript.ffea
```
One final note: for best performance, make sure the number of threads you're running the simulation on divides cleanly into the number of nodes in your rod. For example, if the rod has 14 nodes, you may want to consider running
```bash
export OMP_NUM_THREADS=7
```
## Analysis of rods
`FFEA_rod` comes with its own set of analysis tools, which can be accessed by creating an instance of the `anal_rod` class.

```python
my_analysis = FFEA_rod.anal_rod(my_rod)
```
To compute all of the energies in the rod, run
```python
my_analysis.get_equipartition()
```
This will populate the rod variables `my_analysis.stretch_energy`, `my_analysis.twist_energy` and `my_analysis.bending_energy`. You can also run each test individually `my_analysis.using get_bending_response_mutual()`, `my_analysis.get_twist_amount()` and `my_analysis.get_stretch_energy()`.

You can plot the results of the integration test by running
```python
my_analysis.plot(temp=300)
```
This will write some PDF files into the working directory which are named according to the name of the rod trajectory file used to generate them. If you have a B matrix that is homogenous and isotropic (which you should, because the equiparition test won't work otherwise), it will write another PDF to the working directory comparing the persistence length of the rod trajectory to that of a worm-like chain with a young's modulus computed from the value of the B matrix.

Some other functions that might be helpful:
```python
length = my_analysis.get_absolute_length(0, my_rod.num_elements, 20) # get the absolute length of the 0th frame. The first and second parameters are the first and last nodes, the last parameter is the frame number.

my_analysis.align_to_equil() # align the rod trajectory to the equilibrium structure. Requires the 'icp' module by Clay Flannigan: https://github.com/ClayFlannigan/icp

node_rmsd = my_analysis.get_node_rmsd() # get the per-node rmsd between the current and best fit structures. If you want best-fit rmsd, run align_to_equil first or run with the parameter 'align=True'

time_rmsd = my_analysis.get_time_rmsd() # get the rmsd of all the nodes against time. If the trajectory is already aligned, pass this the parameter 'is_aligned=True'.

eigvals = my_analysis.get_B_eigenvalues() # get the eigenvalues of the B matrix.

my_analysis.thin(100) # thin out the trajectory to have a particular number of frames.
```

In addition to the functions listed above, you can also write your own! All of the data is stored as numpy arrays, so it should be easy enough. You might notice that the FFEA_rod class doesn't store or calculate the values of p, only the values of r. To compute p, just use
```python
my_analysis.p_i = my_analysis.get_p_i(rod.current_r)
my_analysis.equil_p_i = my_analysis.get_p_i(rod.equil_r)
```

Note: by default, the FFEA_rod module will try to compile some math functions into a shared library and import them with ctypes. If it can't, it will fall back to (slower) python implementations. If it has failed, the variable 'FFEA_rod.rod_math_core_status'. will be false.  

This software is still in relatively early development. If you encounter problems, please leave an issue on our [issue tracker](https://bitbucket.org/FFEA/ffea/issues?status=new&status=open).  For more detailed information on the code behind the rod algorithm, please visit the developer documentation [todo: link here].

## Creating rods from all-atom trajectories

You can create an FFEA_rod file from an all-atom trajectory. The method behind this is described in the paper 'KOBRA: A Fluctuating Elastic Rod Model for Slender Biological Macromolecules', which isn't out yet. In summary, this method selects rod parameters such that the local mean square fluctuations in
bending, twisting and extension will match those observed in the MD trajectory. As you might expect, this needs a relatively long trajectory to be accurate. It also needs the trajectory to be composed solely of coiled-coils, so if your trajectory contains globular regions, cut those out!

This tool is accessible in two ways - from the new FFEA python API, and from the ffeatools command-line. For the API, it is accessible under the following namespace:

```python
import ffeatools
rod_analysis, delta_omega, L_i, B = ffeatools.modules.rod.cc_extractor(...)
```

Note: you can't use this module in the old-style (e.g. `import FFEA_rod`) FFEA module imports because they aren't structured to allow importing of modules outside of the FFEA `modules` path.

The ndc_extractor function returns an instance of the rod analysis object, which also contains an FFEA_rod object (`rod_analysis.rod`). It also returns the values used to compute the anisotropic B matrix.

To use this tool from the terminal:
```bash
ffeatools ndc_extractor [args]
```
The parameters are the same for both versions and are accessible under the `ndc_extractor` docstring or with `ffeatools ndc_extractor --help`. They are:

* `prmtop_file` - path to the prmtop (topology) file generated by AMBER.
* `mdcrd_traj_file` - path to the mdcrd trajectory file generated by AMBER.
* `inpcrd_file` - path to the inpcrd file generated by AMBER. Note: this will be used as the equilibrium configuration for the rod, and if it's not actually the equilibrium, all the rod constants will be wrong.
* `rod_out` - path to the output rod file to be written.
* `unroll_rod` - whether to 'unroll' the rod before computing the various constants. This means setting up the rod such that the equilibrium material axes are all pointing the same way. For technical reasons, this can make B easier to compute.
* `get_B` - whether to also compute the anisotropic\inhomogeneous B matrices
* `get_beta` - whether to compute inhomogeneous values for beta, the twist constant.
* `get_kappa` - whether to compute inhomogeneous values for kappa, the stretching constant.
* `target_length` - the number of elements in the rod to be created. This should be high enough to capture all the dynamics of the molecule, but no higher. The default is 15.
* `cluster_size` - the size of the cluster of atoms used to define a node. This 'clustered averaging' is used to avoid averaging out useful dynamics. The default is 10.
* `radius` - the radius, in meters, of the coiled-coil. Most are about 5e-9m.

## Creating a rod-blob coupling

Rods can connect to FFEA blobs in order to form systems containing a mixture of coiled-coils and globular domains, which are extremely common. In order to have a rod-blob coupling, an FFEA simulation must contain at least one rod and at least one blob. The parameters of the coupling are set by a 'coupling' block in the FFEA file, and each block refers to only a single coupling (e.g. a rod with one blob at each end will have two connections).

To start with, create an FFEA file populated with all the rods and blobs that you might want. Then, for each coupling, add a block in the .ffea file that looks like this:

```xml
	<coupling type = blob-to-rod>
        <rod_id = 0>
        <blob_id = 0>
        <rod_node_id = 0>
        <blob_element_id = 0>
        <blob_node_ids = (0,1,2)>
        <node_weighting = (-1,-1,-1)>
        <rotation = (0,0,0)>
        <order = 0>
    </coupling>
```

The parameters should be set as follows:

* `coupling type` - can be either blob-to-rod or rod-to-blob. This is a purely technical distinction, but you can think of all rods as starting at element 0 and going to element n. Going from 0 to n is 'forward' and going from n to 0 is 'backward' along the rod. If the connection is at the 'end' of the rod (element n) this should be set to 'rod-to-blob'. If it at the start of the rod, it should be set to 'blob-to-rod'.
* `rod_id` - the numerical id of the rod to connect. These are indexed by the order that they appear in your .ffea file, starting at zero.
* `blob_id` - the same, but for blobs.
* `rod_node_id` - the id of the node that the interface is at (either 0 or n - you can put the connection in the middle, if you really want, but I wouldn't recommend it).
* `blob_element_id` - the id of the element (the tetrahedron) on the blob that the rod is going to connect to. I would suggest finding the nodes first, then using FFEA_rod.rod_creator.get_connection_info to get the id of the element containing those nodes. Check the function docstring for more information.
* `blob_node_ids` = the ids of the three linear nodes on the element to which the rod will connect. You can easily view these by opening up the FFEA viewer, setting 'add atoms' to 'onto linear nodes', and typing 'label all, resi' into the pymol console.
* `node_weighting` - how to position the connection inside the node. The three values are between 0 and 1 and correspond to the three tetrahedron edges eminating from the internal node (the node which isn't on the face). The resulting position will be the sum of these edges multiplied by these three values. In general, either (-1, -1, -1) or (0.3333333, 0.3333333, 0.3333333) should be enough, that will put the coupling node in the middle of the element
* `rotation` - euler rotation angles relative to the normal of the surface face of the connection element. In radians.
* `order` - this just refers to the order in which the connection is set up. When you create a rod-blob coupling, the equilibrium position of the coupled object becomes relative to the object that it's coupling to. For example, if you create a connection that goes from a blob to a rod, then the starting position of the rod will be relative to the position of the blob. The 'order' is the order that these equilibrium positions are calculated in. If you are creating a chain of objects, the first object in the chain should have order 0, the second should have order 1, etc. Note: you can set this to -1 if you'd rather set up your connections manually, but be wary - the simulation will crash if you start with connections way out of equilibrium.

There are a few other things to consider when setting up rod-blob connections. Make sure that the simulation box is large enough contain the longest axis of your chain. Don't set up a simulation with a connection crossing over a periodic boundary. If possible, use hard boundary conditions.

This feature could be considered to be in early alpha, and its mechanics are rather opaque, but a publication should be out soon. If you have questions or comments, please pay a visit to the FFEA issue tracker. Hypothetical future editions of the software might also make setting up a connection less cumbersome. 
