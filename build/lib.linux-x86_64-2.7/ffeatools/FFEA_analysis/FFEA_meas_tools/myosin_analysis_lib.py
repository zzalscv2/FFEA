# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 15:52:30 2015

@author: py12rw
"""

import FFEA_trajectory, FFEA_pin, FFEA_script
import numpy as np
import sys
import matplotlib.pyplot as plot

def run_sequential_tests(input_ffea_file,
                         head_pin_file,
                         tail_pin_file,
                         head_point_index,
                         tail_point_index,
                         center_point_index,
                         head_line_indices=None,
                         tail_line_indices=None,
                         twist_test=True,
                         dist_test=True,
                         angular_distribution=True,
                         save_plots = True,
                         molecule_name="Myosin 7"):
    """
    This is a library for analysing data from Myosin-7 trajectory files (but it
    could easily be used on any type of molecule with an applicable structure).
    The tests it performs are:
    - End-to-end distance for each frame
    - Angular distribution for each frame
    - Twist (in rads/m) for each frame
    ---
    In:
    input_ffea_file, the path to the .ffea script file associated with your
    desired trajectory.
    head_pin_file: an FFEA_pin file containing a list of indices of nodes in
    the head of the molecule.
    tail_pin_file: the same but in the tail.
    head_point_index: the index of a single node in the center of the head.
    center_point_index: the index of a single node in the center (ideally in
    the bendiest part).
    tail_point_index: the index of a single node in the center of the tail.
    head_line_indices: a list containing the indices of a number of points
    in the head, which form a straight line orthogonal to the long axis of
    the molecule.
    tail_line_indices: the same but for the tail - should be pointing the same
    direction as head_line_indices, though.
    twist_test: boolean, whether to perform a twist test
    dist_test: boolean, whether to perform an end-to-end distance test. Note:
    if this is false, then the script will plot twist angle (rads) instead of
    twist (rads/m).
    angular_distribution: boolean, whether to perform an angular distribution
    test
    save_plots: boolean,  whether to plot the results and save them to PDF
    files in the working directory.
    molecule_name: string, containing the name of the molecule (used in plots).
    Out:
    A dictionary file with named arrays containing the results of each test.
    """
                             
    print("Loading in files...")
    script = FFEA_script.FFEA_script(input_ffea_file)
    trajectory = script.load_trajectory()
    head_pin = FFEA_pin.FFEA_pin(head_pin_file)
    tail_pin = FFEA_pin.FFEA_pin(tail_pin_file)
    head_subblob = trajectory.blob[0][0].define_subblob(head_pin.get_node_indices())
    tail_subblob = trajectory.blob[0][0].define_subblob(tail_pin.get_node_indices())
    
    print("Testing...")    
    
    results = {}
    results["dist_trajectory"] = np.array([])
    results["lever_angle_trajectory"] = np.array([])
    results["twist_angles"] = np.array([])
    results["twist_amount"] = np.array([])
    
    if dist_test:
        results["dist_trajectory"] = trajectory.calc_distance_between_subblobs(0, 0, head_subblob, tail_subblob)
    
    if twist_test:
        results["twist_angles"] = twist_analysis(trajectory, head_pin, tail_pin, head_line_indices, tail_line_indices)

    if angular_distribution:
        results["lever_angle_trajectory"] = trajectory.get_lever_angle_trajectory(head_point_index, center_point_index, tail_point_index)

    if twist_test and dist_test:
        results["twist_amount"] = results["twist_angles"]/results["dist_trajectory"]

    if save_plots:
        filename = input_ffea_file.split(".")[0]
        if "/" in filename:
            filename = filename.split("/")[-1]
        time_axis = get_time_axis(script, trajectory)
        plot_results(results, time_axis, filename, molecule_name)

    return results

def get_time_axis(script, trajectory):
    """
    Gets an array with which to populate the x axis of plots instead of the
    less-useful frame numbers. Assumes a constant timestep.
    In: script (FFEA_script object) and trajectory
    (FFEA_trajectory.FFEA_trajectory) object.
    Out: a 1-d array.
    """
    dt = script.params.dt
    frames_axis = np.arange(0, trajectory.num_frames)
    time_axis = frames_axis*dt
    return time_axis

def plot_results(results, time_axis, filename, molecule_name):
    """
    Plots the results of the end-to-end-distance test, lever angle test and
    twist_amount test, if they exist. Saves pdf plots to the working directory.
    In: results, a dictionary containing a number of named 1-d arrays created
    in run_sequential_tests(), time_axis, the times of each framenumber (also
    created in run_sequential_tests), the string to prefix the filenames with,
    and a string containing the molecule name (shows up on the plots).
    Out: nothin'. Just side effects.
    """
    plot.xlabel("Time (s)")
    
    if results["dist_trajectory"].any():
        plot.xlabel("Time (s)")
        plot.plot(time_axis, results["dist_trajectory"])
        plot.title(molecule_name+" End-to-end Distance")
        plot.ylabel("End-to-end distance (m)")
        plot.savefig(filename+"_endtoend.pdf", format="pdf")
        plot.cla()
        
    if results["lever_angle_trajectory"].any():
        plot.xlabel("Time (s)")
        plot.plot(time_axis, results["lever_angle_trajectory"])
        plot.title(molecule_name+" Lever Angle Evolution")
        plot.ylabel("Lever angle (rads)")
        plot.savefig(filename+"_leverangle.pdf", format="pdf")
        plot.cla()
        
    if results["twist_amount"].any():
        plot.xlabel("Time (s)")
        plot.plot(time_axis, results["twist_amount"])
        plot.title(molecule_name+" Twist Amount")
        plot.ylabel("Twist amount (rads/m)")
        plot.savefig(filename+"_twist_amount.pdf", format="pdf")
        plot.cla()
    elif results["twist_angles"].any():
        plot.xlabel("Time (s)")
        plot.plot(time_axis, results["twist_angles"])
        plot.title(molecule_name+" Twist Angle")
        plot.ylabel("Twist angle (rads)")
        plot.savefig(filename+"_twist_angles.pdf", format="pdf")
        plot.cla()
    
    plot.close()
    return

def twist_analysis(trajectory, head_pin, tail_pin, head_line_indices, tail_line_indices):
    """
    Perform twist analysis on a trajectory.
    For a given trajectory object, this will work out the 'twist' in the
    molecule for each frame. It needs pins containing all the nodes in the head
    and tail, plus a list of the indices defining two lines in the head and
    tail, in the same direction but orthogonal to the long axis of the
    molecule. This can be the .get_node_indices() of a pinfile, if you please,
    but either way, you'll have to get the indices of the nodes manually at
    some point.
    In: trajectory object (type of FFEA_trajectory.FFEA_trajectory), two
    pinfiles (type of FFEA_pin.FFEA_pin) and two lists.
    Out: a 1-d array containing twist angles, in radians.
    """
    #load in data

    if type(head_line_indices) == type(None) or type(tail_line_indices) == type(None):
        raise Exception("Error: no indices supplied for head\tail lines")

    print("Initialising...")

    head = trajectory.blob[0][0].define_subblob(head_pin.get_node_indices())
    tail = trajectory.blob[0][0].define_subblob(tail_pin.get_node_indices())

    print("Getting centroids...")    
    
    head_centroid = trajectory.blob[0][0].get_centroid_trajectory(head)
    tail_centroid = trajectory.blob[0][0].get_centroid_trajectory(tail)
    #optimisation note: loading in the centroids is by far the slowest part
    #of the whole script, but if I calculate the centroids beforehand and then
    #try and pass them to this function, I get a weird bug/side effect and
    #the centroid disappears. There's probably an easy way to make it work,
    #though.    
    
    twist_angles = np.zeros(len(head_centroid))

    print("Projecting and calculating angles...")    
    
    for frame_num in range(len(head_centroid)):
        
        head_point = head_centroid[frame_num]
        tail_point = tail_centroid[frame_num]
        
        head_line_transformed = np.zeros([len(head_line_indices), 3])
        tail_line_transformed = np.zeros([len(tail_line_indices), 3])
        
        # transfooorm!

        for i in range(len(head_line_indices)):
            head_line_transformed[i] = transform_point(trajectory.blob[0][0].frame[frame_num].pos[head_line_indices[i]], head_point, tail_point)
        for i in range(len(tail_line_indices)):
            tail_line_transformed[i] = transform_point(trajectory.blob[0][0].frame[frame_num].pos[tail_line_indices[i]], head_point, tail_point)
            #optimisation note: this bit used to use two different planes instead
            #of one - you could make it faster by rewriting transform_point to
            #trasform two points in one plane.

        head_vector = get_vector_from_average_points(head_line_transformed)
        tail_vector = get_vector_from_average_points(tail_line_transformed)
        
        angle = get_angle_between_vectors(head_vector, tail_vector)
        
        twist_angles[frame_num] = angle
         
    print("Done calculating twist angles!")    
    
    return twist_angles

##################
#                #
#  VECTOR STUFF  #
#                #
##################

def transform_point(point_to_transform, head_centroid, tail_centroid):
    """
    This is a bit of a wrapper function that uses the stuff below to project
    an arbitrary point into an arbitrary plane. The plane is defined using two
    points on the normal to that plane (the head and tail centroids).
    In: point_to_transform, head_centroid_tail_centroid: 1-d numpy arrays of
    3 elements each.
    Out: a point that has been transformed into this plane, also as a 1-d array.
    """
    plane_normal = get_vector_from_points(head_centroid, tail_centroid)
    plane_center = get_midpoint(head_centroid, tail_centroid)
    t = get_t_from_plane_normal_and_position(plane_normal, plane_center, point_to_transform)
    transformed_point = get_new_point_coords(plane_normal, point_to_transform, t)
    return transformed_point

def get_midpoint(point1, point2):
    return (point1+point2)/2.0
    
def get_vector_from_points(point1, point2):
    """
    Returns a direction vector from two points.
    In: two points, as 1-d numpy arrays with 3 elements
    Out: the direction vector, also as a 1-d numpy array with 3 elements.
    """
    return (point1 - point2)
    
def get_t_from_plane_normal_and_position(plane_normal, plane_position, point):
    """
    Calculates the length of a vector t that goes between the surface of a
    plane and an arbitrary point, in the direction of the normal to that plane.
    In: plane normal, plane position and point, three 1-d numpy arrays with 3
    elements each.
    Out: t, a scalar.
    """
    x, y, z = point[0], point[1], point[2]
    a, b, c = plane_normal[0], plane_normal[1], plane_normal[2]
    d, e, f = plane_position[0], plane_position[1], plane_position[2]
    t = (a*d - a*x + b*e  - b*y + c*f - c*z)/(a**2 + b**2 + c**2)
    return t
    
def get_new_point_coords(plane_normal, point, t):
    """
    Calculates the position of a point in an arbitrary plane, given t, the
    distance of the point from the plane in the direction of the normal
    defining that plane.
    In: plane_normal and point, both 1-d numpy arrays with 3 elements, and t,
    which is just a scalar (should be a float)
    Out: a 1-d numpy array with the co-ordinates of the point in the plane.
    """
    return point + plane_normal*t
    
def get_vector_from_average_points(points_array): # a 2-d array with cols for x y and z
    """
    Gets an average from a 2-d array that lists a series of points in 3D space.
    Returns the average x, y, and z co-ordinate (or whatever floats your boat)
    in a 1-d numpy array.
    """
    return np.array([np.sum(points_array[:,0]), np.sum(points_array[:,1]), np.sum(points_array[:,2])]  )
    
def get_angle_between_vectors(vec1, vec2): #it's sad that these are for 3d points only they should be n-dimensional really
    """
    Gets the angle between two vectors in radians.
    In: vec1, vec2, both numpy array objects with 3 elements.
    Out: an angle in radians.
    """
    dot_product = np.sum(vec1*vec2)
    mag1 = np.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
    mag2 = np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
    return np.arccos(dot_product/(mag1*mag2))