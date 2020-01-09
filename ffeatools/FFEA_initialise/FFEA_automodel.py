# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

import subprocess
import argparse as _argparse
import FFEA_node as _FFEA_node, FFEA_topology as _FFEA_topology, FFEA_surface as _FFEA_surface
import sys
try:
    import __builtin__
except ImportError:
    import builtins as __builtin__
    print("\n---\nWarning: you are running ffea using Python 3. Not every FFEA script has been updated to be forward-compatible yet - if you experience errors, please try downgrading to Python 2!")

#import math
import os

# Set up argparse
parser = _argparse.ArgumentParser(description="Create an FFEA model and script for an electron density map.")
parser.add_argument("--params", action="store", dest='params', type=str, default=["1.5e3", "1e-3", "1e-3", "5.5e8", "2.2e9"], nargs="+", help="Material parameters (density, shear viscosity, bulk viscosity, shear modulus, bulk modulus). The defaults are the values for water (1.5e3 1e-3 1e-3 5.5e8 2.2e9).")
parser.add_argument("--granularity", action="store", dest='granularity', type=str, default="10", help="Granularity of the model in Angstroms. Set this to the size of the smallest feature you want to still be resolvable. Default value is 12.") #args.o
parser.add_argument("--cull", action="store", dest='cull', type=str, default="10000", help="Remove small floaters and holes from the model. In cubic Angstroms. Default value is 10000.") #args.o
parser.add_argument("iso", action="store", type=str, default="1.15", help="Set ISO level for EM density map.") #args.o
parser.add_argument("--threshold", action="store", dest='threshold', type=float, default=20.0, help="Area of smallest permissible face. The default is 25 square angstroms. Be careful of lowering this value, as it could make your simulations unstable!") #args.o
parser.add_argument("file", action="store", type=str, help="Input electron density map (.map) file.")

def get_exitcode_stdout_stderr(args):
    """
    Execute the external command and get its exitcode, stdout and stderr.
    Will throw an OSError if that thing does not exist.
    """
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode
    proc.wait()
    return exitcode, out, err
    
def check_error(code, out, err):
    """
    Check for an error code returned by a subprocess. If such an error is
    returned, show it.
    """
    if code != 0 or err:
        print(out)
        if err:
            raise RuntimeError(err)
        else:
            raise RuntimeError("\n---\nCould not create model.")    

def get_area_heron(thing):
    """
    Use heron's formula to get the area of a triangle with vertices specified
    in 3D space. The trinagle is given as a list of lists.
    """
    r1, r2, r3 = thing[0], thing[1], thing[2]
    abs_a = ( (r1[0]-r2[0])**2+ (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2 )**0.5
    abs_b = ( (r2[0]-r3[0])**2+ (r2[1]-r3[1])**2 + (r2[2]-r3[2])**2 )**0.5
    abs_c = ( (r3[0]-r1[0])**2+ (r3[1]-r1[1])**2 + (r3[2]-r1[2])**2 )**0.5
    s = (abs_a+abs_b+abs_c)/2
    return ( s*(s-abs_a)*(s-abs_b)*(s-abs_c) )**0.5
    
def check_stl(filename, threshold):
    """
    Parse the (poorly-formatted) STL file and check that no face has an area
    below a certain threshold (in angstroms squared).
    """
    # Load stl
    vertex_data = []
    with open(filename) as fp:
        for line in fp:
            if line.split(" ")[0] == "		vertex":
                line.replace("  "," ")
                line_new = line.split(" ")[1:]
                while '' in line_new:
                    line_new.remove('')
                for i in range(len(line_new)):
                    line_new[i] = float(line_new[i])
                vertex_data.append( line_new )

    for i in range(len(vertex_data)/3):
        A = vertex_data[i*3]
        B = vertex_data[(i*3)+1]
        C = vertex_data[(i*3)+2]
        area = get_area_heron([A, B, C])
        if area<threshold:
            return False
        
    return True
    
def get_E(shear_mod, bulk_mod):
    """
    Get the young's modulus (E) given the shear and bulk moduli.
    """
    return (9*bulk_mod*shear_mod)/( (3*bulk_mod)+shear_mod)
    
def validate_args(args):
    """
    Perform some simple sanity checks on the argparse args list.
    """
    
    if sys.version_info[0] > 2:
        print("\n---\nWarning: you are running ffea using Python 3. Not every FFEA script has been updated to be forward-compatible yet - if you experience errors, please try downgrading to Python 2!")

    if os.path.exists(args.file) == False:
        raise IOError("\n---\nCan't find a file at "+args.file+".")
    
    if float(args.granularity)<3:
        print("\n---\nWarning: your granularity is set very fine. This may generate very small elements. It will be automatically increased if the elements are found to be too small.")
    if float(args.granularity)>30:
        print("\n---\nWarning: your granularity is set very coarse. After the model has been created, please check to make sure it looks right, or consider decreasing the value.")
    
    if get_E(float(args.params[3]), float(args.params[4])) > 15000000000:
        print("\n---\nWarning: your material is set to be very stiff. If it's not supposed to be, please change your material parameters, and\or use a long timestep!")

    if get_E(float(args.params[3]), float(args.params[4])) < 100000000:
        print("\n---\nWarning: your material is set to be very flexible. If it's not supposed to be, please change your material parameters, and\or use a short timestep!")

    try:
        tetgentest_code, tetgentest_out, tetgentest_err = get_exitcode_stdout_stderr(["tetgen"])
    except OSError:
        raise OSError("\n---\ntetgen is not installed, but is needed for mesh generation. If you don't have a module or package, you can download the tetgen source from http://wias-berlin.de/software/tetgen/. If you have a tetgen binary, please make sure it has been appended to your system $PATH!")

    try:
        ffeatest_code, ffeatest_out, ffeatest_err = get_exitcode_stdout_stderr(["ffea", "-h"])
    except OSError:
        raise OSError("\n---\nFFEA is not installed. Please use the FFEA installer script SCRIPT NAME\PATH HERE before creating a model.")

    try:
        toolstest_code, toolstest_out, toolstest_err = get_exitcode_stdout_stderr(["ffeatools", "--help"])
    except OSError:
        raise OSError("\n---\nFFEAtools are not installed, or are not on the system path. If you see this error, either you've run this script before installing ffeatools, or your install is broken. If you're in the latter camp, please leave an issue on our issue tracker!")

    if args.file.split(".")[-1] == "pdb":
        raise IOError("\n---\nPlease use an electron density map file (.map), not a pdb. If you only have a PDB, you can convert it using the following command:\n\n  ffeatools pdbtoemmap [input pdb filename] [output map filename] [nx] [ny] [nz] [atomic radius]\n\nWhere nx, ny and nz represent the voxel resolution of your desired map, and atomic radius is the radius of a sphere centered on each atom.")

    return
    
def check_obj(obj_out):
    """
    Check that there are a minimum amount of stuff written into the .obj. If not, then the ISO level is wrong.
    (obj_out is a string containing the output of emmaptosurf)
    """
    for line in obj_out.split("\n"):
        linesplit = line.split(" ")
        if linesplit[0] == "Wrote":
            if int(linesplit[1]) < 1000 or int(linesplit[4]) < 1000:
                return False
            else:
                return True
            
def check_danger_elements(vol):
    """
    From the volume data file, this function tells the user if there
    are any extremely long, thin elements in the simulation. If so, it returns
    True, if not, false.
    """
    top = _FFEA_topology.FFEA_topology()
    top.load_vol(vol)
    node = _FFEA_node.FFEA_node()
    node.load_vol(vol)
    
    eindex = 0
    dindex = []
    for e in top.element:
        if e.get_smallest_lengthscale(node) < 0.5: #assuming .vol in nm
            print e.get_smallest_lengthscale(node)
            dindex.append(eindex)
        eindex += 1
    if len(dindex)>0:
        return False
    return True
        
def automodel(args):
    """
    This script will check initialisation parameters and do boring stuff for you.
    It's not intended for experenced FFEA users, rather new ones.
    Features:
        - Automatically runs all of the FFEA init scripts with sensible defaults
        - Checks that the right tools are installed and prompts the user if not
        - Checks that the ISO level setting produced an actual mesh
        - Checks that the material parameters are specified in the right units and aren't simulation-destroying
        - Checks the granularity, and auto-adjusts it if small elements are created
        - Prompts the user with options when they're done initialising
    This function takes an argparse args list as its parameters, see the top of the file for more info on those.
    """
    
    print("Validating input...")
    validate_args(args)

    name = args.file.split(".")[0]
    
    print("Cleaning up electron density map...")
    surf_code, surf_out, surf_err = get_exitcode_stdout_stderr(["ffeatools", "emmaptosurf", "-map", args.file, "-out", name+"_processed.map", "-format", "map", "-level", args.iso, "-cull_floaters", args.cull, "-fill_cavities", args.cull])
    check_error(surf_code, surf_out, surf_err)  
    if "Couldn't extract data" in surf_out:
        raise IOError("Unable to read EM map file.")
    
    print("Converting .map to .obj...")
    obj_code, obj_out, obj_err = get_exitcode_stdout_stderr(["ffeatools", "emmaptosurf", "-map", name+"_processed.map", "-out", name+".obj", "-format", "obj", "-level", args.iso])
    check_error(obj_code, obj_out, obj_err)
    obj_OK = check_obj(obj_out)
    if obj_OK == False:
        raise IOError("\n---Incorrect scale for the ISO level.")
    
    attempts = 1
    stl_correct = False
    no_danger_elements = False
    current_granularity = float(args.granularity)
    
    while stl_correct == False or no_danger_elements == False:
        print("Coarse graining (attempt "+str(attempts)+" at granularity "+str(current_granularity)+"A)...")
        grain_code, grain_out, grain_err = get_exitcode_stdout_stderr(["ffeatools", "surftocgsurf", name+".obj", name+".stl", str(current_granularity), "y", "y"])
        check_error(grain_code, grain_out, grain_err)
        if check_stl(name+".stl", args.threshold):
            stl_correct = True
        
        print("Filling in tetrahedra with tetgen...")
        tet_code, tet_out, tet_err = get_exitcode_stdout_stderr(["tetgen", "-Y", name+".stl"])
        check_error(tet_code, tet_out, tet_err)
    
        print("Converting to .vol format...")
        vol_code, vol_out, vol_err = get_exitcode_stdout_stderr(["ffeatools", "tettonet", "-i", name+".1.ele", name+".1.face", name+".1.node", "-o", name+".vol"])
        check_error(vol_code, vol_out, vol_err)
        
        
        no_danger_elements = check_danger_elements(name+".vol")
        
        attempts+=1
        current_granularity+=1.
    
    try:
        import random
        if random.random()>0.97:
            print("Cloning Monsters and GROwing shrimps...")
    except:
        pass
    

    
    print("Creating FFEA script file...")
    
    try:
        import ffeatools # this one is in the API!
        sys.stdout = open(os.devnull, "w") # supress printing (yes, this is cross-platform!)
        ffeatools.FFEA_initialise.FFEA_convert_from_volume.convert_from_volumetric_mesh(name+".vol", stokes_radius=None, cull=0.0, density=float(args.params[0]), shear_visc=float(args.params[1]), bulk_visc=float(args.params[2]), shear_mod=float(args.params[3]), bulk_mod=float(args.params[4]), dielectric=1.0, make_script=True, outfname=None)
        sys.stdout = sys.__stdout__
    except ImportError: # but installing the API is optional so...
        ffea_code, ffea_out, ffea_err = get_exitcode_stdout_stderr(["ffeatools", "voltoffea", "--mesh", name+".vol", "--density", args.params[0], "--shear-visc", args.params[1], "--bulk-visc", args.params[2], "--shear-mod", args.params[3], "--bulk-mod", args.params[4], "--dielec", "1.0", "--make-script"])
        check_error(ffea_code, ffea_out, ffea_err)
    
    print("---")
    print("Finished!")
    print("")
    print("To tweak your simulation parameters (e.g. number of steps, timestep, kT, solvent stokes viscosity), open up your FFEA script file, "+name+".ffea.")
    print("")
    print("To take a look at your model, open up your script file in the FFEA PyMOL plugin (you can find the plugin in the FFEA source folder, under /ffeatools/analysis/)")
    print("")
    print("Finally, to run your simulation, type 'ffea "+name+".ffea' (without quotes)")
    sys.stdout = open(os.devnull, "w") # supress printing
    return

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    # This will run the contents of the script, unless it is being run from an IDE or from inside the FFEA API.
    args = parser.parse_args()
    # Get args and build objects
    automodel(args)
