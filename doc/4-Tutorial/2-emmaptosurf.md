EDM to Surface Profile {#emmaptosurftut}
=============================

Once we have a piece of volumetric data such as an em density map, we can create a surface profile of it. However, we first need to make sure that the density map is completely filled, and doesn't contain random offshoots that we don't want (artifacts of the experimental data for example). Our EM density map, as viewed from the top down, contains density corresponding to the GroEl molecule being filled when, atomistically, we know that this is not the case. We can account for this by firstly slightly raising the recommended countour level, to separate this mass:

![The electron density structure of EMDB ID:EMD-5043 (countour level 1.15).](emd5043_2.png "GroEL Electron Microscopy Structure")

`ffeatools` contains a program which is a modified version of the flood fill algorithm used in graphics software, to fill internal cavities and get rid of separated section:

	ffeatools emmaptosurf -map [INPUT .map fname] -out [OUTPUT fname] -format [OUTPUT format (map/obj/surf)] -level [isolevel] -cull_floaters [Voxel threshold] 
			          -fill_cavities [Voxel threshold] -pdb [OPTIONAL INPUT original .pdb fname]

The relevent options are cull_floaters, which gets rid of separated mass smaller than the given threshold, -fill_cavities, which fills in cavities within the structure smaller than the given threshold, and -level, which is the isolevel we were working at within Chimera. As you can see, we can output immediately to another electron density map to check our progress, and so if we use the following parameters:
	
	ffeatools emmaptosurf -map emd_5043.map -out emd_5043_processed.map -format map -level 1.15 -cull_floaters 10000 
			          -fill_cavities 10000

We get the following structure:

![The electron density structure of EMDB ID:EMD-5043 post-processing (countour level 1.15).](emd5043_processed_both.png "GroEL Electron Microscopy Structure")

We can now use exactly the same program to generate a surface profile from our processed EM map, in either .stl or .obj format:

	ffeatools emmaptosurf -map emd_5043_processed.map -out emd_5043.obj -format obj -level 1.15

Visualising using Netgen 6.0 (http://sourceforge.net/projects/netgen-mesher/):

![STL Surface profile of EMDB ID:EMD-5043.](emd5043_stl.jpg "GroEL Surface Profile")
	
