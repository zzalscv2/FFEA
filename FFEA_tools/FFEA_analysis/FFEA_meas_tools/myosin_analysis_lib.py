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
                         center_pin_file,
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
    center_pin = FFEA_pin.FFEA_pin(center_pin_file)
    head_subblob = trajectory.blob[0][0].define_subblob(head_pin.get_node_indices())
    tail_subblob = trajectory.blob[0][0].define_subblob(tail_pin.get_node_indices())
    
    print("Testing...")    
    
    results = {}
    results["dist_trajectory"] = np.array([])
    results["lever_angle_trajectory"] = np.array([])
    results["head_angles"] = np.array([])
    
    if dist_test:
        results["dist_trajectory"] = trajectory.calc_distance_between_subblobs(0, 0, head_subblob, tail_subblob)
    
    if twist_test:
        results["head_angles"], results["tail_angles"] = twist_analysis(trajectory, head_pin, tail_pin, center_pin, head_line_indices, tail_line_indices)

    if angular_distribution:
        results["lever_angle_trajectory"] = trajectory.get_lever_angle_trajectory(head_point_index, center_point_index, tail_point_index)

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
        plot.plot(time_axis, results["dist_trajectory"])
        plot.title(molecule_name+" End-to-end Distance")
        plot.ylabel("End-to-end distance (m)")
        plot.savefig(filename+"_endtoend.pdf", format="pdf")
        plot.cla()
        
    if results["lever_angle_trajectory"].any():
        plot.plot(time_axis, results["lever_angle_trajectory"])
        plot.title(molecule_name+" Lever Angle Evolution")
        plot.ylabel("Lever angle (rads)")
        plot.savefig(filename+"_leverangle.pdf", format="pdf")
        plot.cla()
        
    if results["head_angles"].any():
        plot.plot(time_axis, results["head_angles"], label="Head")
        plot.plot(time_axis, results["tail_angles"], label="Tail")
        plot.title(molecule_name+" Twist Angle")
        plot.ylabel("Angle (rad)")
        plot.legend()
        plot.savefig(filename+"_twist.pdf", format="pdf")
        plot.cla()
        
    return

def twist_analysis(trajectory, head_pin, tail_pin, center_pin, head_line_indices, tail_line_indices):

    #load in data

    if type(head_line_indices) == type(None) or type(tail_line_indices) == type(None):
        raise Exception("Error: no indices supplied for head\tail lines")

    print("Initialising...")

    head = trajectory.blob[0][0].define_subblob(head_pin.get_node_indices())
    tail = trajectory.blob[0][0].define_subblob(tail_pin.get_node_indices())
    center = trajectory.blob[0][0].define_subblob(center_pin.get_node_indices())

    print("Getting centroids...")    
    
    head_centroid = trajectory.blob[0][0].get_centroid_trajectory(head)
    tail_centroid = trajectory.blob[0][0].get_centroid_trajectory(tail)
    center_centroid = trajectory.blob[0][0].get_centroid_trajectory(center)
    
    head_angles = np.zeros(len(head_centroid))
    tail_angles = np.zeros(len(head_centroid))

    print("Projecting and calculating angles...")    
    
    for frame_num in range(len(head_centroid)):
        
        head_point = head_centroid[frame_num]
        tail_point = tail_centroid[frame_num]
        center_point = center_centroid[frame_num]
        
        head_line_transformed = np.zeros([len(head_line_indices), 3])
        tail_line_transformed = np.zeros([len(tail_line_indices), 3])
        
        # transfooorm!

        for i in range(len(head_line_indices)):
            head_line_transformed[i] = transform_point(trajectory.blob[0][0].frame[frame_num].pos[head_line_indices[i]], head_point, center_point)
        for i in range(len(tail_line_indices)):
            tail_line_transformed[i] = transform_point(trajectory.blob[0][0].frame[frame_num].pos[tail_line_indices[i]], center_point, tail_point)
        
        #project into 2d        
        
        head_line_transformed_projection = head_line_transformed[:,0:2] # get only the 1st 2 dimensions
        tail_line_transformed_projection = tail_line_transformed[:,0:2] # change this if that's wrong
        
        frame_head_gradient = get_gradient_from_points(head_line_transformed_projection)
        frame_tail_gradient = get_gradient_from_points(tail_line_transformed_projection)
        
        if frame_num==0:
            head_angles[0] = 0
            first_head_gradient = frame_head_gradient
            first_tail_gradient = frame_tail_gradient
        else:
            head_angles[frame_num] = get_angle_between_two_gradients(first_head_gradient, frame_head_gradient)
            tail_angles[frame_num] = get_angle_between_two_gradients(first_tail_gradient, frame_tail_gradient)
    
    print("Done!")    
    
    return head_angles, tail_angles

################
#              #
#  MISCELLANY  #
#              #
################

def get_angle_between_two_gradients(m1, m2): #in radians
    return np.arctan( (m2-m1)/(1 + (m2*m1))  )
    
def get_gradient_from_points(points): # a 2-d array wherein each point has its own row
    A = np.rot90(np.array([points[:,0], np.ones(len(points[:,0]))]), 3) #why numpy wants it in this format i have no idea
    y = points[:,1]
    m, c = np.linalg.lstsq(A, y)[0]
    return m

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