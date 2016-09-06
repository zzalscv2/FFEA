import sys, os
import FFEA_script
from FFEA_universe import *
import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Convert .vol files to the standard gamut of FFEA files (.ffea, .top, .surf, .mat, etc)")
parser.add_argument("--mesh", action="store", help="Input mesh file in .vol format")
parser.add_argument("--density", action="store", type=float, default=1.5e3, help="Density of material (kg/m^3)")
parser.add_argument("--shear_visc", action="store", type=float, default=1e-3, help="Shear viscosity of material (Pa.s)")
parser.add_argument("--bulk_visc", action="store", type=float, default=1e-3, help="Bulk viscosity of material (Pa.s)")
parser.add_argument("--shear_mod", action="store", type=float, default=370e6,help="Shear modulus of material (Pa)")
parser.add_argument("--bulk_mod", action="store", type=float, default=111e7, help="Shear viscosity of material (Pa)")
parser.add_argument("--stokes_radius", action="store", help="Stokes radius (for hydrodynamics)")
parser.add_argument("--dielec", action="store", type=float, default=1.0, help="Dielectric constant (unitless)")
parser.add_argument("--make_script", action="store_true", help="Whether to generate a .ffea script file")
parser.add_argument("--cull", action="store", type=float, help="Cull all elements smaller than a certain volume")

def convert_from_volumetric_mesh(mesh, stokes_radius=None, cull=[False, 0.0], density=1.5e3, shear_visc=1e-3, bulk_visc=1e-3, shear_mod=370e6, bulk_mod=111e7, dielectric=1.0, make_script=False):
    """
    This script converts the .vol file from the initialisation routines into the necessary file formats for an FFEA simulation
    This means get linear element from vol and make 2nd order, store new faces, tets etc, build vdw, pin, bsites, stokes, move blob to centroid
    
    # Node - 2nd order
    # Top - 2nd order
    # Surf - 2nd order
    # Mat - 1st order
    # Pin - N/A or 2nd order (linked to nodes)
    # VdW - 2nd Order
    # lj - N/A
    # Stokes - 2nd order
    # Bsites - N/A or 2nd order (linked to surf)
    # beads - N/A
    """    
    # Get args
    volfname = mesh
    
    # Test args
    if volfname == None:
        raise IOError("You must specify a .vol filename with '--mesh'. Run this script with the -h flag to get help.")
    
    basename = os.path.splitext(os.path.abspath(volfname))[0]
        
    # Get initial stuff from .vol file!
    node = FFEA_node.FFEA_node(volfname)
    top = FFEA_topology.FFEA_topology(volfname)
    surf = FFEA_surface.FFEA_surface(volfname)
    
    # Let each surface face know which element it is connected to (if this is slow, just load surface from topology instead)
    if surf.get_element_indices(top) == -1:
        surf = top.extract_surface()
    
    # Now, build necessary things that only have linear properties
    mat = FFEA_material.FFEA_material()
    mat.build(top.num_elements, d=density, sv=shear_visc, bv=bulk_visc, sm=shear_mod, bm=bulk_mod, di=dielectric)
    
    # Convert stuff to 2nd order
    top.increase_order(node=node, surf=surf)
    
    # Find out what is interior and what is surface and reorder stuff
    top.calculateInterior(surf=surf)
    node.calculateInterior(top=top, surf=surf)
    
    # Check the normals in the and surface files only (this is all that is necessary for FFEA)
    surf.check_normals(node, top)
    
    # Cull small elements
    if cull:
        top.cull_interior(cull, node, surf=surf)
    
    # Everything should be set! Make default files for the other stuff
    pin = FFEA_pin.FFEA_pin()
    vdw = FFEA_vdw.FFEA_vdw()
    vdw.default(surf.num_faces)
    lj = FFEA_lj.FFEA_lj()
    lj.default()
    
    # User may set a stokes radius. If they don't, make a default one
    if stokes_radius == None:
    
        # All drag is currently local, so get largest length scale and treat as sphere. Very bad approximation for long molecules
        dims = node.calculate_dimensions()
        stokes_radius = max(dims) / len(top.get_linear_nodes())
    
    stokes = FFEA_stokes.FFEA_stokes()
    stokes.default(node.num_nodes, top, stokes_radius)
    
    # Script will be defaulted
    if make_script:
        script = FFEA_script.FFEA_script()
        script.default(basename)
    
    # Now, print them all out!
    node.write_to_file(basename + ".node")
    top.write_to_file(basename + ".top")
    surf.write_to_file(basename + ".surf")
    mat.write_to_file(basename + ".mat")
    vdw.write_to_file(basename + ".vdw")
    lj.write_to_file(basename + ".lj")
    pin.write_to_file(basename + ".pin")
    stokes.write_to_file(basename + ".stokes")
    
    if make_script:
        script.write_to_file(basename + ".ffea")

if sys.stdin.isatty():
    args = parser.parse_args()
    convert_from_volumetric_mesh(args.mesh, args.stokes_radius, args.cull, args.density, args.shear_visc, args.bulk_visc, args.shear_mod, args.bulk_mod, args.dielec, args.make_script)
