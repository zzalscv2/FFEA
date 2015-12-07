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
#                         center_pin_file,
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
    dt = script.params.dt
    frames_axis = np.arange(0, trajectory.num_frames)
    time_axis = frames_axis*dt
    return time_axis

def plot_results(results, time_axis, filename, molecule_name):
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
        
    return

def twist_analysis(trajectory, head_pin, tail_pin, head_line_indices, tail_line_indices):

    #load in data

    if type(head_line_indices) == type(None) or type(tail_line_indices) == type(None):
        raise Exception("Error: no indices supplied for head\tail lines")

    print("Initialising...")

    head = trajectory.blob[0][0].define_subblob(head_pin.get_node_indices())
    tail = trajectory.blob[0][0].define_subblob(tail_pin.get_node_indices())

    print("Getting centroids...")    
    
    head_centroid = trajectory.blob[0][0].get_centroid_trajectory(head)
    tail_centroid = trajectory.blob[0][0].get_centroid_trajectory(tail)
    
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

def transform_point(point_to_transform, end_centroid, center_centroid):
    plane_normal = get_vector_from_points(end_centroid, center_centroid)
    t = get_t_from_plane_normal_and_position(plane_normal, center_centroid, point_to_transform)
    transformed_point = get_new_point_coords(plane_normal, point_to_transform, t)
    return transformed_point

def get_midpoint(point1, point2):
    return (point1+point2)/2.0
    
def get_vector_from_points(point1, point2):
    return (point1 - point2)
    
def get_t_from_plane_normal_and_position(plane_normal, plane_position, point):
    x, y, z = point[0], point[1], point[2]
    a, b, c = plane_normal[0], plane_normal[1], plane_normal[2]
    d, e, f = plane_position[0], plane_position[1], plane_position[2]
    t = (a*d - a*x + b*e  - b*y + c*f - c*z)/(a**2 + b**2 + c**2)
    return t
    
def get_new_point_coords(plane_normal, point, t):
    return point + plane_normal*t
    
def get_vector_from_average_points(points_array): # a 2-d array with cols for x y and z
    return np.array([np.sum(points_array[:,0]), np.sum(points_array[:,1]), np.sum(points_array[:,2])]  )
    
def get_angle_between_vectors(vec1, vec2): #it's sad that these are for 3d points only they should be n-dimensional really
    dot_product = np.sum(vec1*vec2)
    mag1 = np.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
    mag2 = np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
    return np.arccos(dot_product/(mag1*mag2))