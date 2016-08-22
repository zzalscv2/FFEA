Surface Profile to Volumetric Mesh {#surftovoltut}
=============================

Once we have our structure, it is a `simple' matter of filling the surface with tetrahedra. This can be done either with NETGEN or TETGEN.

NETGEN Meshing
=============

Meshing with NETGEN allows for a direct visualisation of the structure as it is being meshed. However, it is a little tricky to install and the meshing options are not at all clear. The process is as follows:

  * Open NETGEN
  * Load the surface. Click File->Load Geometry... and select a .stl file
  * Choose the meshing parameters. Click Mesh->Mesh Options. Vary the mesh granularity for an approximate set of options, or click the Mesh Size tab for more control
  * Fill the structure with tetrahedra. Select 'Generate Mesh' and wait until it finishes
  * Save the volume file. Click File->Save Mesh... and save the file as .vol (we'll call ours emd_5043_8ang.vol)

![Loading a .stl file using NETGEN](netgenloadstlprintscreen.png "NETGEN - Loading an STL Surface")
<BR>
![Changing the meshing options in NETGEN](netgenoptionsprintscreen.png "NETGEN - Meshing parameters interface")
<BR>
![Constructing a .vol file using NETGEN](netgengenmeshprintscreen.png "NETGEN - Building a Volume")
<BR>
![Saving a .vol file using NETGEN](netgensavemeshprintscreen.png "NETGEN - Saving a Volume")
<BR>
The .vol file contains all of the structural information we need to build an FFEA system.
![Coasened volume profiles of EMDB ID:5043 built using NETGEN. From the top, moving clockwise: 1A, 3A, 5A, 8A, 10A](netgencoarsening.png "NETGEN - Coarsening Process")

TETGEN Meshing
=============

Meshing with TETGEN is by default a simple command line interface. However, in this developers opinion, the optional flags allow for a much better control over the final meshed structure. Therefore, we also provide a conversion script to transform the TETGEN output files into a NETGEN .vol format for visualisation purposes (this was much easier than I'm making it out to be...). The process is as follows:


	tetgen -Y emd_5043_8ang.stl

As simple as that. The -Y flags keeps the surface structure exactly as it was, ensuring elements are optimised without compromising their lower size limit as much as possible. We now have a series of files; emd_5043_10ang.1.* that contain the volumetric data. Let's convert this into the NETGEN .vol format for completeness:

	python /path/to/FFEAtools/FFEA_initialise/Volume_tools/convert_tet_to_net.py -i emd_5043_8ang.1.ele -i emd_5043_8ang.1.face -i emd_5043_8ang.1.node -o emd_5043_8ang.vol

![Coasened surface profiles of EMDB ID:5043 built using TETGEN. From the top, moving clockwise: 1A, 3A, 5A, 8A, 10A](tetgencoarsening.png "TETGEN - Coarsening Process")

Now we have all of the relevent structural data, all that is left is to convert it into the formats required for an FFEA simulation! Whether you prefer NETGEN or TETGEN is completely up to the you, so long as we are left with a .vol file at the end, FFEA will run perfectly well!
