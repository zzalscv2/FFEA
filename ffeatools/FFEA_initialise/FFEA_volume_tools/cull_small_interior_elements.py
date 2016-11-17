import sys, os
from FFEA_universe import *
import FFEA_script
import numpy as np
import argparse as _argparse
import __builtin__

parser = _argparse.ArgumentParser(description="Cull small interior elements")
parser.add_argument("input_file", action="store", help="Input file (.vol or .ffea)")
parser.add_argument("smallest_length", action="store", help="Lower length limit)")

def cull_small_interior_elements(input_file, smallest_length, load_file=False):

    # Get args
    base, ext = os.path.splitext(input_file)
    print ext
    
    node = None
    top = None
    surf = None
    
    if load_file:
    
        if ext == ".vol":
        
            try:
                node = FFEA_node.FFEA_node(sys.argv[1])
                top = FFEA_topology.FFEA_topology(sys.argv[1])
                smallest_length = float(sys.argv[3])
            except:
                sys.exit("Could not load the input structure from '.vol' file.")
        
        elif ext == ".ffea":
            try:
                script = FFEA_script.FFEA_script(sys.argv[1])
                node = script.load_node(0)
                top = script.load_topology(0)
                surf = script.load_surface(0)
            except:
                sys.exit("Could not load the input structure from '.ffea' file.")
                
    else:
        node = input_file.load_node(0)
        top = input_file.load_topology(0)
        surf = input_file.load_surface(0)
        
    # Firstly, check which elements are linear
    
    # Do we already know?
    for i in range(top.num_elements):
    
        #if top.element[i].interior == None:
        
        # Let's assume that we don't know, so calculate from surface (fast) or from topology (slow)
        calc_from_topology = False
        if surf != None:
    
            # Test if surface knows about it's topology
            for j in range(surf.num_faces):
                if surf.face[j].elindex == -1:
                    calc_from_topology = True
                    break
    
        if calc_from_topology:
    
            # Slow. Calculate from topology
            for j in range(top.num_elements):
                sys.stdout.write("\rTesting whether element %d of %d is interior (First element loop will be slow, sorry :()" % (j, top.num_elements))
                sys.stdout.flush()
                top.isElementInterior(j)
        else:
            # Fast. Read from surface
            # Set all to interior first
            for j in range(top.num_elements):
                top.element[j].interior = True
                
            for j in range(surf.num_faces):
                top.element[surf.face[j].elindex].interior = False
        break
    
    #for i in range(top.num_elements):
        #if top.element[i].interior:
        #    print i
    sys.exit()
    # Let's do some coarsening
    current_smallest_length = 0.0
    
    # While limit not satisfied
    print "Coarsening...\n"
    while current_smallest_length < smallest_length:
        
        current_smallest_length = float("inf")
        for i in range(top.num_elements):
            sys.stdout.write("\rShould we coarsen element %d of %d?)" % (i, top.num_elements))
            sys.stdout.flush()
    
            # Only consider interior elements (surface elements would be well complex, and change the overall shape)
            if top.isElementInterior(i):
                
                # Don't just get any element for length < smallest_length. Get current smallest element, for stability reasons
                length = np.power(6 * top.element[i].get_volume(node), (1.0/3.0))
                if length < current_smallest_length:
                    print "\r%d is interior with a length scale of %f                        " % (i, length)
                    current_smallest_length = length
                    delete_index = i
    
        print delete_index, current_smallest_length
        
        # Delete this element by creating a node at it's centre, reconnecting all nodes to that point. 4 nodes should be deleted, 1 element deleted and 1 node added
    
        # Create node and add to list
        new_node_pos = top.element[i].calc_centroid(node)
        node.add_node(new_node_pos, nodetype = 1)
    #return whatever the result was
    
if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    # Get args and build objects

    cull_small_interior_elements(args.input_file, args.smallest_length, load_file=True)
    #save result to file here