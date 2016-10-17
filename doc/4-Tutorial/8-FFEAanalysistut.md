FFEA Analysis {#FFEAanalysistut}
=============================

## Setting up

In previous parts of this tutorial, we have employed several python scripts by invoking 'ffeatools' at the terminal. These scripts all make use of a set of core FFEA python modules. Each module corresponds to a different FFEA data file - the trajectory (``.out``), pin files (``.pin``), the stokes file (``.stokes``), et cetera. These are all unified and linked by the FFEA script file (``.ffea``).

In addition to a few basic analysis tools that can be run from the terminal, these core Python modules are provided as part of a Python package, which allows them to be imported into a Python interpreter. To install this package, open the FFEA source folder (containing ``setup.py``) in the terminal, and type

	python setup.py install

This will install the FFEA tools into your Python site-packages folder. For this tutorial, you may want to consider installing the [Anaconda Python distribution](https://www.continuum.io/downloads), which comes with a set of common scientific Python packages, and allows for packages to be installed and removed easily. An interactive Python shell (such as IPython) or a Python IDE with an object inspector and auto-complete is recommended. Anaconda comes with one called Spyder, which can be launched using the terminal

	spyder

In our interactive session we first
```python
import ffeatools
```
The FFEA tools can now be accessed in the ffeatools namespace. To start, we create a new script object.

```python
>>> our_script = ffeatools.modules.script("/path/to/script.ffea")
```
## The FFEA Trajectory

We can now load in the contents of any of the files that the script is linked to. We will start by loading the trajectory.

```python
>>> our_trajectory = our_script.load_trajectory()
Loading FFEA trajectory file...
Frames read = 1, Frames skipped = 0
Frames read = 2, Frames skipped = 0
Frames read = 3, Frames skipped = 0
...
Frames read = 99, Frames skipped = 0
Frames read = 100, Frames skipped = 0
Frames read = 101, Frames skipped = 0
done! Successfully read 101 frame/s from '/path/to/trajectory.out'.
```
We can now access all of the data inside the trajectory. For example,

```python
>>> our_trajectory.num_blobs
1
```

The trajectory object has a hierarchical structure:

* Trajectory object
  * Meta information
  * Methods
  * Blobs list
    * Blob 0
      * Conformations list
        * Conformation 0
          * Centroid
          * Frames list
            * Frame
              * Node position array
        * Conformation 1 (etc...)
        * Centroid
    * Blob 1 (maybe!)
    * Blob 2 (maybe!)
    * etc...

Each FFEA trajectory is composed of one or more blobs (the example in the previous tutorial, GroEL, is just one blob). These blobs are then separated info conformations. The number of nodes in each conformation is the same, but the  position of these nodes is different. Each conformation contains a list of all the frames in that trajectory, and inside each of these frames is a two-dimensional array containing the positions of all of the nodes in the trajectory.

>A brief aside: This documentation currently does not cover switching the conformation of a blob during a simulation, as the feature is still under development. For that reason, the conformation index should always be 0. If a molecule does have multiple conformations in the same simulation, then both conformations will contain all the frames, and the node positions of frames in inactive conformations is '``None``'.

The data in the trajectory can be accessed like this:

```python
>>> our_trajectory.blob[blob_index][subblob_index].frame[frame_num].pos[node_index]
```

For example, the zeroth node of the zeroth frame of the zeroth subblob of the zeroth blob is:

```python
>>> our_trajectory.blob[0][0].frame[0].pos[0]
array([  2.65931600e-09,  -2.47283800e-09,  -1.88632000e-09])
```

A one-dimensional array containing the positions of those nodes.

Each blob contains a 'centroid'. The centroid trajectory is the trajectory of the average position of every node in the blob. This can be extremely useful, as individual nodes often fluctuate so much that analysing data from them is extremely difficult. To get the centroid of a blob,

```python
>>> our_centroid_trajectory = our_trajectory.blob[0][0].get_centroid_trajectory()
```

Will return an array, containing arrays of the sub-blob indices and the trajectory itself. For example, the first frame of the centroid trajectory is

```python
>>> our_centroid_trajectory[1][0]
array([  1.46518608e-17,   5.21608637e-18,  -4.92614644e-18])
```

## Sub-blobs and Pinfiles
Blobs can also contain sub-blobs. By default, the zeroth sub-blob of a blob contains the indices of all the nodes in that blob. The user can create additional sub-blobs which only contain a subset of these indices, and use that sub-blob to get a trajectory or a centroid.

To create a new sub-blob, we specify the indices of the points that we wish to use. The easiest way to do this is to create a pin. A pin is a file or object containing a list of indices. They have many functions in FFEA, but they were primarily designed to specify groups of nodes that have a different behaviour (for example, nodes that are anchored in place). To create a pin, we spin up an instance of the pin module:

```python
>>> our_pin = ffeatools.modules.pin()
	File '' not found.
```

This will create an empty pin (supplying a filename can load in a pre-existing pin). We can fill it with the indices of the nodes we wish to pin manually...

```python
>>> for index in [1, 2, 3, 4, 5]:
>>>   our_pin.add_pinned_node(index)
```

Or we can pin large groups of nodes at once, using the radial pin function. To use this function, we first specify a blob and conformation to use:

```python
>>> pin_frame = our_trajectory.blob[0][0].frame[0]
```

Then we supply this to the pin_radially method, along with the index of the central node and the radius of the pin in meters.

```python
>>> our_pin.pin_radially(pin_frame, 1, 0.000000001)
```

We can also request the indices of the nodes in the pin.

```python
>>> our_pin.index
[1, 11, 14, 321, 375, 470, 600, 615, 691, 702, 810, 1067, 1095, 1232, 1656]
```

We can now create a subblob using this pin by running

```python
>>> our_trajectory.blob[0][0].set_subblob(our_pin)
```

And get its trajectory with

```python
>>> our_trajectory.blob[0][0].get_centroid_trajectory(1)
```

Where 1 is the index of our subblob. If we were to create another subblob, it would be given the index 2, et cetera. Remember that the zeroth index is reserved for a subblob containing the indices of every node.

To save a pin for later, use
```python
>>> our_pin.write_to_file("our_pin.pin")
```

## Material and Vdw Files

In performing analysis on trajectories, you may want to extract the data from other FFEA files. Most of these files are structured in a very similar way. For example, the .mat files contain information about the material parameters of each element, and can be acccessed like this:

```python
>>> our_material = our_script.load_material # create an instance
>>> our_material(0,0).element[0] # load the material parameters from the 0th element
array([  1.50000000e+03,   1.00000000e-03,   1.00000000e-03,
         5.50000000e+08,   2.20000000e+09,   1.00000000e+00])
```
The elements in the array, from first to last, are: the density of the material, in \f$kg/m^3\f$, the shear viscosity, in \f$Pa \cdot s\f$, the bulk viscosity, in \f$Pa \cdot s\f$, the shear modulus, in \f$Pa\f$, the bulk modulus, in \f$Pa\f$, and the dielectric constant, which is unitless.

Similarly, to get the Van Der Waals type for each face:

```python
>>> our_vdw = our_script.load_vdw(0,0) # (blob 0, conformation 0)
Loading FFEA vdw file...
>>> our_vdw.index[0] # get the vdw type of the 0th face
-1
```

## Topographies and Surfaces

At this point, you may be wondering what the use for the Van Der Waals type for a particular face is if we don't know which nodes and which elements are associated with that face. As a reminder:

* The .node file tells us the positions, in 3-D space, of each node (a point)
* The .top file tells us the connectivity of these nodes, how they are arranged into elements (tetrahedrons)
* The .surf file tells us the connectivity of the surface nodes: how they are arranged into faces (triangles)

> Aside: If we have a toplogy, why even have a surface at all? When we generate these files, we start with the surface, and we fill in the volume with tetgen. But, at certain points during the initialisation process, before we have a topology, we still have to deal with the surface - hence, the 'surface' module. Later, we end up with a topology, and the topology contains all of the information about the surface, and more. You can easily go from a topology to a surface - you just have to work out which faces are being shared between two elements, and which ones are on their own. This can be a bit slow, so we keep the surface file around anyway, in case we need it for something.

Loading the node file results in something that ought to be familiar:

```python
>>> our_node = our_script.load_node(0,0) #(0th blob, 0th conformation)
Loading FFEA node file...
>>> our_node.pos[0]
array([ 108.833827,   51.687247,   25.75118 ])
```
Finally, we can load the topology in like so:

```python
>>> our_topology = our_script.load_topology(0,0)
Loading FFEA topology file...
```

The topology is made up of elements, and the elements can be accessed like this:
```python
>>> our_topology.element[0].n # find the 0th element
[859, 887, 32, 4, 1540, 1254, 1474, 1403, 1316, 466]
```

> Another aside: why are there 10 values? Isn't this element supposed to comprise a tetrahedron? The first four elements are indeed the elements that make up the tetrahedron, but the next six actually make up points on the same tetrahedron.

> If we think about just one face of the tetrahedron, that face is described by three points, connected by straight (linear) lines. But some calculationas in FFEA (such as electrostatics) can make use of second-order elements - the lines connecting the sides of the triangle are no longer linear, they can bend inward or outward. This second set of values are at the midpoint between two 'first-order' nodes, and thus describe the second-order behaviour of the element. In most simualtions (and most analyses) they can be safely ignored.

In the FFEA_topography module, we can also calculate the volume of a given element. As the topology module only contains information about the connectivity of the nodes, we need to supply our node object in order to retrieve the volume:

```python
>>> our_topology.element[0].calc_volume(our_node)
208.98603491219549
```

And we can compare two elements to see if they share an edge:

```python
>>> our_topology.element[0].sharesanedge(our_topology.element[1])
False
```
