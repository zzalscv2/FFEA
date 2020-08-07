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

import argparse as _argparse
import FFEA_node as _FFEA_node, FFEA_topology as _FFEA_topology, FFEA_surface as _FFEA_surface
import sys
try:
    import __builtin__
except:
    import builtins as __builtin__

from os.path import splitext
from os.path import exists

# Set up argparse
parser = _argparse.ArgumentParser(description="Convert tetgen files (.ele, .face, .node) to a .vol file (FFEA, netgen)")
parser.add_argument("-i", action="store", nargs="+", help="Input files (.ele, .face, .node).")
parser.add_argument("-o", action="store", help="Output .vol file.")

def convert_tet_to_net(input_fnames, output_fname):
    """
    Converts tetgen-compatible files to a netgen-compatible file (which is also
    compatible with FFEA).
    In:
    input fnames : a list of filenames or paths to files for .ele (topology),
    .face (surface) and .node files in the tetgen format.
    output_fname : the name of the desired output file e.g. "output.vol"
    Out:
    nothing in python, but it dumps a file in the working directory.
    """

    # Read files
    top = None
    surf = None
    node = None
    
    if type(input_fnames) == type(None) or (len(input_fnames) != 3 and len(input_fnames) != 4):
        print("Expected 3 input files! Use --help for help.")
        raise IOError

    for f in input_fnames:
        base = f.rsplit(".", 1)[0]
        ext = f.rsplit(".", 1)[1]
        
        if ext == "ele":
            top = _FFEA_topology.FFEA_topology(f)
    
            # Additionally
            if output_fname == "":
                output_fname = base + ".vol"
    
        elif ext == "node":
            node = _FFEA_node.FFEA_node(f)
        elif ext == "face":
            surf = _FFEA_surface.FFEA_surface(f)
    
    if output_fname == None:
        print("\tOutput filename not found. Creating default output filename...")
        output_fname = splitext(input_fnames[0])[0] + ".vol"
        if exists(output_fname):
             print("\tDefault output filename already exists. Please supply your own filename (which will be overwritten if it exists")
             raise IOError

    # Now, write them all to file
    fout = open(output_fname, "w")
    fout.write("mesh3d\ndimension\n3\ngeomtype\n11\n\n")
    fout.close()
    
    surf.write_to_file(output_fname)
    top.write_to_file(output_fname)
    node.write_to_file(output_fname)
    
    fout = open(output_fname, "a")
    fout.write("endmesh\n")
    fout.close()

# Run argparse if we detect a tty
if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
        convert_tet_to_net(args.i, args.o)
    except IOError:
        parser.print_help()
