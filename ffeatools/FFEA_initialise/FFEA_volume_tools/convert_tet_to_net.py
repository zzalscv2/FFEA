import argparse as _argparse
import FFEA_node as _FFEA_node, FFEA_topology as _FFEA_topology, FFEA_surface as _FFEA_surface
import sys

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
    
    if type(input_fnames) == type(None) or len(input_fnames) != 3:
        raise IOError("Expected 3 input files! Use --help for help.")
    
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
if sys.stdin.isatty():
    args = parser.parse_args()
    convert_tet_to_net(args.i, args.o)