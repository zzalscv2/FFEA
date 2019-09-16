# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:53:47 2017

@author: rob
"""
import icp
import FFEA_pdb
import numpy as np
import FFEA_script
import sys
import argparse as _argparse
from scipy.linalg import expm
import copy
import random
import os

try:
    import __builtin__
except ImportError:
    import builtins as __builtin__
    print("\n---\nWarning: you are running ffea using Python 3. Not every FFEA script has been updated to be forward-compatible yet - if you experience errors, please try downgrading to Python 2!")
    
parser = _argparse.ArgumentParser(description="Rotate an FFEA model onto a PDB file using the iterative closest point method.")
parser.add_argument("--bindex", action="store", dest='bindex', type=int, default=0, help="If you have multiple blobs, give the index of the blob that you want to align here.") #args.o
parser.add_argument("--cindex", action="store", dest='cindex', type=int, default=0, help="If you have multiple conformations, give the conformation that you want to align here.") #args.o
parser.add_argument("script", action="store", type=str, default="1.15", help="The script that corresponds to the FFEA structure to align.") #args.o
parser.add_argument("--iterations", action="store", dest='iterations', type=int, default=2000, help="Max number of iterations of the alignment algorithm before stopping. The default is 2000.") #args.o
parser.add_argument("--tolerance", action="store", dest='tolerance', type=float, default=0.00001, help="The minimum amount of improvement in RMSD needed to keep the alignment algorithm running. The default is 0.00001, which already a bit excessive. ") #args.o
parser.add_argument("--candidates", dest='candidates', action="store", type=int, default=100, help="Number of different starting locations for the alignment algorithm. The default is 50, but increase this value if it looks like the algorithm gets stuck before it finds the true minimum rmsd.") #args.o
parser.add_argument("pdb", action="store", type=str, help="PDB file to align the FFEA structure to.")
parser.add_argument("--node", dest='node', action="store_true", default=False, help="Add this flag to load and align a node file (for alignment before you run your simulation).") #args.o
parser.add_argument("--traj", dest='traj', action="store_true", default=False, help="Add this flag to load and align a trajectory file (for if you've already run the simulation).") #args.o
parser.add_argument("--no_save", action="store_true", dest='no_save', default=False, help="Do not modify the FFEA files. This will only print the rotation matrix.")

def rot_euler(v, xyz):
    ''' Rotate vector v (or array of vectors) by the euler angles xyz '''
    # https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    for theta, axis in zip(xyz, np.eye(3)):
        v = np.dot(np.array(v), expm(np.cross(np.eye(3), axis*-theta)))
    return v

def align_centroid(node_object, pdb_object):
    """
    Extract the positional daat from the two objects. Create a centroid for the
    pdb object (the entire object doesn't have one, only the different chains)
    and return thevector that puts the node object onto the PDB.
    """
    pdb_centroid = np.array([0.,0.,0.])
    for chain in pdb_object.chain:
        pdb_centroid += chain.frame[0].calc_centroid()
    pdb_centroid/=len(pdb_object.chain)
    node_centroid = node_object.calc_centroid()
    diff =  pdb_centroid - node_centroid
    node_object.translate(diff)
    return diff
    
def create_atom_array(pdb_object):
    """
    Given a PDB object, return an array with the position of every atom inside
    (again, there is normally just one for eacch chain.)
    """
    atom_no = 0
    chain_no = 0
    num_atoms = pdb_object.num_atoms[0]
    atom_array = np.zeros([num_atoms,3])
    while atom_no < len(atom_array):
        try:
            atom_array[atom_no] = pdb_object.chain[chain_no].frame[0].pos[atom_no]
            atom_no += 1
        except IndexError:
            chain_no += 1
    return atom_array

def apply_transformation_4x4(pos_array, T):
    """
    Apply a 4x4 transformation matrix to our 2-D array of 3-D points. Returns
    the newly translated array.
    """
    translated_array = np.zeros([ len(pos_array), 3 ])
    for i in range(len(pos_array)):
        node_1 = [pos_array[i][0], pos_array[i][1], pos_array[i][2], 1]
        translated_array[i] = np.dot(T, node_1)[:3]
    return translated_array

def align_traj_frame(traj_frame, node_centroid, T):
    """
    Makes up for the discrepancy between the scale\centroids of the node and
    trajectory file, so that the same transformation as the node file can be
    applied. This modifies the trah_frame object directly, there is no
    return value.
    """
    scale_factor = 10**10
    traj_frame_centroid = traj_frame.calc_centroid()*scale_factor
    diff = traj_frame_centroid - node_centroid
    traj_frame.pos = traj_frame.pos*scale_factor - diff
    traj_frame.pos = apply_transformation_4x4(traj_frame.pos, T)
    traj_frame.pos = traj_frame.pos/scale_factor
    
def fit_from_candidates(node_array, pdb_array, max_iterations, req_tolerance, num_candidates):
    """
    This program uses an ICP (iterative closest point) algorithm to get a
    transformation matrix that will align the FFEA structure to a PDB. However,
    this algorithm tends to get stuck in local minima. This function will
    randomly rotate the FFEA structure and run the algorithm, and output
    the rotation and the following translation for the result with the lowest
    RMSD.
    """
    candidates= {}
    for candidate in range(num_candidates):
        print_progress_bar("Finding optimal alignment", num_candidates, candidate)
        temp_node_array = copy.copy(node_array)
        XYZ = [random.random()*2*np.pi, random.random()*2*np.pi, random.random()*2*np.pi]
        rot_euler(temp_node_array, XYZ)
        T, distances = icp.icp(temp_node_array, pdb_array, max_iterations=max_iterations, tolerance=req_tolerance)
        candidates[np.average(distances)] = [T, XYZ]
    return min(candidates.items())

def print_progress_bar(start_text, num_iterations, num_done):
    """
    Badly print an ascii progress bar.
    """
    try:
        rows, cols = os.popen('stty size', 'r').read().split()
        cols = int(cols)
    except ValueError:
        cols = 80, 80
        
    ratio_done = float(num_done)/float(num_iterations)
    str_done = str(num_done)+"/"+str(num_iterations)+" "
    textstr = start_text+" ["
    bar_length= cols-len(textstr.decode('utf-8')+"]")-len(str_done.decode('utf-8'))
    blank_chars = "▒"*int((bar_length*(1-ratio_done)))
    filled_chars = "█"*int((bar_length*(ratio_done)))
    spacer = ""
    str_to_write = "\r"+textstr+filled_chars+blank_chars+"] "+str_done
    
    if len(str_to_write.decode('utf-8')) > cols:
        blank_char_num_new = int((bar_length*(1-ratio_done)))-(len(str_to_write.decode('utf-8'))-cols)
        blank_chars = "▒"*blank_char_num_new
        str_to_write = "\r"+textstr+filled_chars+blank_chars+"] "+str_done
        
    if len(str_to_write.decode('utf-8')) < cols:
        spacer = " "*(cols-len(str_to_write.decode('utf-8')))
        
    if len(str_to_write.decode('utf-8')) != cols:
        return

    sys.stdout.write(str_to_write+spacer)

def main(script_file, pdb_file, num_iterations=2000, req_tolerance=0.00001, no_save=False, bindex=0, conf=0, node=False, traj=False, num_candidates=100):
    """
    Align an FFEA script file to a PDB file.
    Parameters:
        - script_file - path to the .ffea script file
        - pdb_file - path to the PDB file
        - num_iterations - max iterations of the ICP algorithm
        - req_tolerance - convergence criteria for the ICP algorithm
        - no_save - whether to modify the files or just output the matrices
        - bindex - index of the blob to align
        - cindex - index of the conformation to align
        - node\traj- which object to apply the transformation to
    Returns
        - Translation vector, euler rotation, transformation matrix, RMSD
    """
    if node==False and traj==False:
        print("Nothing to do! Please add the parameter --node or --traj, depending on which file you want to align.")
        return
    
    print("Loading stuff...")
    
    if node:
        
        save_stdout = sys.stdout
        sys.stdout = open('trash', 'w')
    
        script = FFEA_script.FFEA_script(script_file)
        node = script.load_node(bindex, cindex=conf)
        pdb = FFEA_pdb.FFEA_pdb(pdb_file)
        
        sys.stdout = save_stdout
        
        diff = align_centroid(node, pdb)
        node.translate(diff)
        
        node_array = node.pos
        pdb_array = create_atom_array(pdb)
        #print("Finding optimal alignment...")
    
        results = fit_from_candidates(node_array, pdb_array, num_iterations, req_tolerance, num_candidates)
        print(" ")
        rmsd = results[0]
        T = results[1][0]
        XYZ = results[1][1]
        
        print("Optimal alignment found.")
        print("Translation: "+str(diff))
        print("Rotation: "+str(XYZ))
        print("Transformation matrix: \n"+str(T))
        print("New RMSD: +"+str(rmsd))
    
        if no_save:
            print("Done! (didn't save anything)")
            return diff, XYZ, T, rmsd
        
        print("Applying rotation...")
        rot_euler(node.pos, XYZ)
        
        print("Applying transformation...")
        node.pos = apply_transformation_4x4(node.pos, T)
        
        print("Saving node file...")
        node.write_to_file(script.blob[bindex].conformation[conf].nodes)
        
        print("Done!")
        
        return
        
    if traj:
        
        save_stdout = sys.stdout
        sys.stdout = open('trash', 'w')
        
        script = FFEA_script.FFEA_script(script_file)
        traj = script.load_trajectory()        
        pdb = FFEA_pdb.FFEA_pdb(pdb_file)

        sys.stdout = save_stdout
        
        #print("Finding optimal alignment...")
        
        for frame_no in range(len(traj.blob[bindex][conf].frame)):
            traj.scale(10**10, frame_no)
        
        diff = align_centroid(traj.blob[bindex][conf].frame[0], pdb)
        
        for frame_no in range(len(traj.blob[bindex][conf].frame)):
            if frame_no>0:
                traj.blob[bindex][conf].frame[frame_no].translate(diff)
                
        pdb_array = create_atom_array(pdb)
        node_array = traj.blob[bindex][conf].frame[0].pos
        
        results = fit_from_candidates(node_array, pdb_array, num_iterations, req_tolerance, num_candidates)
        print(" ")
        rmsd = results[0]
        T = results[1][0]
        XYZ = results[1][1]
        
        print("Optimal alignment found.")
        print("Translation (angstroms): "+str(diff))
        print("Rotation (radians): "+str(XYZ))
        print("Transformation matrix: \n"+str(T))
        print("New RMSD (angstroms): +"+str(rmsd))
        
        if no_save:
            print("Done! (didn't save anything)")
            return diff, XYZ, T, rmsd
        
        print("Applying rotation...")
        
        for frame in traj.blob[bindex][conf].frame:
            rot_euler(frame.pos, XYZ)
        
        print("Applying transformation...")
        
        for frame in traj.blob[bindex][conf].frame:
            frame.pos = apply_transformation_4x4(frame.pos, T)
            
        print("Saving trajectory...")
        
        save_stdout = sys.stdout
        sys.stdout = open('trash', 'w')
        
        for frame_no in range(len(traj.blob[bindex][conf].frame)):
            traj.scale(10**-10, frame_no)
        
        traj.write_to_file(script.params.trajectory_out_fname)

        sys.stdout = save_stdout
        
        print("Done!")
        print("Note: if you load this trajectory into the viewer, the first frame will load from the node file instead, so it will look like it isn't rotated. PLEASE CHECK THE SECOND FRAME!")
        
        return
        
if __name__ == "__main__" and hasattr(__builtin__, 'FFEA_API_mode') == False and sys.stdin.isatty():
    args = parser.parse_args()
    main(args.script, args.pdb, num_iterations=args.iterations, req_tolerance=args.tolerance, no_save=args.no_save, bindex=args.bindex, conf=args.cindex, node=args.node, traj=args.traj, num_candidates=args.candidates)
