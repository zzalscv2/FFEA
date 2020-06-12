# -*- coding: utf-8 -*-
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

import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import FFEA_rod
import copy
import argparse

parser = argparse.ArgumentParser(description="Create an FFEA_rod structure and extract its parameters from an all-atom simulation.")
parser.add_argument("prmtop_file", action="store", type=str, help="Input .prmtop (AMBER topology) file.")
parser.add_argument("mdcrd_file", action="store", type=str, help="Input .mdcrd (AMBER trajectory) file.")
parser.add_argument("inpcrd_file", action="store", type=str, help="Input .inpcrd (AMBER co-ordinate) file used for equilibrium configuration.")
parser.add_argument("chain1_end", action="store", type=int, help="the index of the atom at which the first chain (the first coil) in the coiled-coil ends (int).")
parser.add_argument("rod_out", action="store", type=str, help="Path to output file to write (.rodtraj type).")
parser.add_argument("--trj_format", action="store", dest='trj_format', type=str, help="format of the trajectory file (string). Uses MDAnalysis formats.")
parser.add_argument("--inp_format", action="store", dest='inp_format', type=str, help="format of the equilibrium structure file (string). Uses MDAnalysis formats.")
parser.add_argument("--unroll_rod", action="store_true", dest='unroll_rod', help="Whether to unroll the rod (force the equilibrium twist angle to be zero). Note: this will not alter the energies, only the relative values.") #args.o
parser.add_argument("--get_B", action="store_true", dest='get_B', help="Whether to compute the values of the B matrix (bending stiffness) for each node") #args.o
parser.add_argument("--get_beta", action="store_true", dest='get_inhomogeneous_beta', help="Whether to compute the values of beta, the twisting stiffness for each node") #args.o
parser.add_argument("--get_kappa", action="store_true", dest='get_inhomogeneous_kappa', help="Whether to compute the values of kappa, the stretching stiffness for each element") #args.o
parser.add_argument("--target_length", action="store", dest='target_length', type=int, default="15", help="Target number of nodes in the rod (default 15).") #args.o
parser.add_argument("--cluster_size", action="store", dest='cluster_size', type=int, default="10", help="Number of atoms which are averaged to create each rod node (default 10).") #args.o
parser.add_argument("--radius", action="store", dest='radius', type=float, default="5e-9", help="Radius of the rod in m (default 5e-9).") #args.o

def try_target_steps(starting_step, init_length, margin):
    
    max_step = int(starting_step+(starting_step*margin*0.5))
    min_step = int(starting_step-(starting_step*margin*0.5))
    
    if min_step < 0:
        min_step = 1
    
    curr_step = min_step + 0
    good_steps = []
    while curr_step < max_step:
        #print("Checking step size "+str(curr_step)+"/"+str(max_step))
        if init_length%curr_step == 0:
            good_steps.append(curr_step)
            curr_step += 1
        else:
            curr_step += 1
    return good_steps

def determine_simplification(arr3d, target_length, margin=0.5):
    """
    For a given array and a target length, create a simplified array which is
    that length. This does a similar thing to simplify(), but it will attempt
    find a step size automatically, and will adjust the target length within
    the margin if the length of the array doesn't divide cleanly by the target
    length.
    Params:
        arr3d - 2d array representing 3d nodes
        target_length - target length of the first axis of the new array
        margin: factor by which the target length can be different
    """
    if len(arr3d[0])%target_length == 0:
        step = len(arr3d[0])/target_length
        return simplify(arr3d, step)
    
    arr3d_length = len(arr3d[0])
    starting_step= int(arr3d_length/target_length)
    
    while True:
        good_steps = try_target_steps(starting_step, arr3d_length, margin)
        if len(good_steps) > 0:
            break
        arr3d_length -= 1

    curr_score = 99e999
    curr_step = None
    for step in good_steps:
        new_score = abs(step - starting_step)
        if new_score < curr_score:
            curr_score = new_score
            curr_step = step
            
    arr3d = arr3d[:,:arr3d_length]
    
    return simplify(arr3d, curr_step)

def simplify(arr3d, step):
    """
    Simplify a 2-d array representing a set of nodes by averaging over clusters
    of sequential nodes.
    Step = number of nodes to average over
    arr3d = rod array, e.g. rod.equil_r
    Returns: simplified array
    """
    if len(arr3d[0])%step != 0:
        raise ValueError("Rod length not divisible by step.")
    shape = np.shape(arr3d)
    shape = np.insert(shape, 2, step)
    shape[1]/=step
    shaped_arr3d = np.reshape(arr3d, shape)
    averaged_arr3d = np.average(shaped_arr3d, axis=2)
    if np.shape(averaged_arr3d)[0] == np.shape(arr3d)[0] and np.shape(averaged_arr3d)[2] == np.shape(arr3d)[2] and np.shape(averaged_arr3d)[1] == np.shape(arr3d)[1]/step:
        return averaged_arr3d
    else:
        raise Exception("Something went wrong. It's Rob Welch's fault. Honestly, this script exists only for his benefit, I don't know why you're using it.")

#def get_anisotropic_twist(analysis, k_BT):
#    """
#    Get the value of the twisting constant, beta, per element.
#    Params:
#        analysis, an instance of FFEA_rod.anal_rod with a trajectory loaded.
#        k_BT, bolztmann's constant in SI units multiplied by the temperature
#        of the system.
#    Returns:
#        anisotropic_twist, an array with the following dimensions: element,
#        average B.
#    """
#    try:
#        analysis.p_i
#    except AttributeError:
#        analysis.p_i = analysis.get_p_i(analysis.rod.current_r)
#        analysis.equil_p_i = analysis.get_p_i(analysis.rod.equil_r)
#    
#    delta_theta_arr = np.zeros( [len(analysis.p_i), len(analysis.p_i[0])-1] )
#    Li = np.zeros( [len(analysis.p_i), len(analysis.p_i[0])-1] )
#    
#    for frame in range(len(analysis.p_i)):
#        for element in range(len(analysis.p_i[frame])-1):
#                    
#            mim1_transported = np.squeeze(np.asarray(analysis.py_rod_math.parallel_transport(analysis.py_rod_math.normalize(analysis.rod.current_m[frame][element]), analysis.py_rod_math.normalize(analysis.p_i[frame][element]), analysis.py_rod_math.normalize(analysis.p_i[frame][element+1] ))))
#            delta_theta = analysis.py_rod_math.get_signed_angle( mim1_transported, analysis.py_rod_math.normalize(analysis.rod.current_m[frame][element+1]), analysis.py_rod_math.normalize(analysis.p_i[frame][element]))
#            equil_mim1_transported = np.squeeze(np.asarray(analysis.py_rod_math.parallel_transport(analysis.py_rod_math.normalize(analysis.rod.equil_m[frame][element]), analysis.py_rod_math.normalize(analysis.equil_p_i[frame][element]), analysis.py_rod_math.normalize(analysis.equil_p_i[frame][element+1] ))))
#            equil_delta_theta = np.squeeze(np.asarray(analysis.py_rod_math.get_signed_angle( equil_mim1_transported, analysis.py_rod_math.normalize(analysis.rod.equil_m[frame][element+1]), analysis.py_rod_math.normalize(analysis.equil_p_i[frame][element]))))
#                
#            delta_theta_arr[frame][element] = delta_theta - equil_delta_theta
#                
#            Li[frame][element] = analysis.py_rod_math.get_length(analysis.equil_p_i[frame][element]) + analysis.py_rod_math.get_length(analysis.equil_p_i[frame][element+1]) 
#                #beta[frame][element] = (Li * k_BT) / ( 2*np.power(delta_theta-equil_delta_theta, 2))
#                #if beta[frame][element] > 1e-28:
#                #    beta[frame][element] = 1e-28
#            
#            if frame == 300 and element == 4:
#                print("mim1_transported: "+str(mim1_transported))
#                print("mim1: "+str(analysis.py_rod_math.normalize(analysis.rod.current_m[frame][element])))
#                print("pi: "+str(analysis.py_rod_math.normalize(analysis.p_i[frame][element+1] )))
#                print("mi: "+str(analysis.py_rod_math.normalize(analysis.rod.current_m[frame][element+1])))
#                print("pim1: "+str(analysis.py_rod_math.normalize(analysis.p_i[frame][element])))
#                print("equil_delta_theta: "+str(np.squeeze(np.asarray(equil_delta_theta))))
#                print("delta_theta: "+str(np.squeeze(np.asarray(delta_theta))))
#                print("Li: "+str(Li[frame][element]))
#                print("B: "+str( (Li[frame][element] * k_BT)/ ( 2*np.power( np.squeeze(np.asarray(delta_theta_arr[frame][element])), 2)) ))
#
#    anisotropic_twist = np.zeros(len(analysis.p_i[0])-1)
#    twist_energy = np.zeros( [len(analysis.p_i), len(analysis.p_i[0])-1] )
#    #for element_no in range(len(anisotropic_twist)): 
#        #anisotropic_twist[element_no] = (np.average(Li[:,element_no]) * k_BT)/ ( 2*np.power(np.average(delta_theta_arr[:,element_no]), 2))
#        
#    B = np.zeros( [len(analysis.p_i), len(analysis.p_i[0])-1] )
#    for frame in range(len(analysis.p_i)):
#        for element in range(len(analysis.p_i[frame])-1):
#            B[frame][element] = (Li[frame][element] * k_BT) / ( np.power( np.squeeze(np.asarray(delta_theta_arr[frame][element])), 2))
#            
#            
#            twist_energy[frame][element] = B[frame][element]/(2*Li[frame][element]) * np.power( np.mod( np.squeeze(np.asarray(delta_theta_arr[frame][element])) + np.pi, 2*np.pi) - np.pi, 2)
#            #twist_energy = beta/(2*Li) * np.power( np.mod(delta_theta-equil_delta_theta + np.pi, 2*np.pi) - np.pi, 2)
#            
#    anisotropic_twist = np.average( np.ma.masked_invalid(B), axis=0)
#
#    return anisotropic_twist, B, twist_energy, delta_theta_arr#, delta_theta_arr, Li
#

def get_inhomogeneous_param(analysis, parm_index, kbT=300*1.38064852e-23):
    # parm_index: 0 = stretch, 1 = twist
    
    params_backup = copy.copy(analysis.rod.material_params)
    
    analysis.rod.material_params[:,:,parm_index] = 1.0
    
    if parm_index == 1:
        analysis.get_twist_amount()
        energies_unit = analysis.twist_energy
    elif parm_index == 0:
        analysis.get_stretch_energy()
        energies_unit = analysis.stretch_energy
    else:
        raise IndexError("Supply the index for a parameter (either 1 (twist) or 0 (stretch))")
    
    #params = 0.5*kbT/(np.average(energies_unit[1:], axis=0))
    
    params, error = get_avg_and_mean_error(energies_unit[1:], kbT*0.5, per_element = True)
    
    analysis.rod.material_params = params_backup
    
    return params, error
    

def get_delta_omega(analysis, fast=False):
    """
    This function computes the value of $\Delta \omega$, the material curvature,
    in the same way that analysis.get_bending_response_mutual() from FFEA_rod does.
    Unfortunately, get_bending_rensponse_mutual shoud really have been a pure
    function and not a method of analysis that doesn't return anything, but
    it's too late for that, I accept my folly. Are you reading this because
    you're interested in how it was implemented, or do you actually have to
    work on this file? Look, what can I say, I made a mistake. It's too late
    for that now. We have to just work with what we've got here.
    Params: analysis, an instance of FFEA_rod.anal_rod
    Returns delta_omega and L_i. two 2-d arrays, specifying these values for
    each frame and node.
    """
    def get_weights(a, b):
        a_length = analysis.py_rod_math.get_length(a)
        b_length = analysis.py_rod_math.get_length(b)
        weight1 = a_length/(a_length+b_length)
        return weight1

    def get_mutual_frame_inverse(pim1, pi, weights=None):
        if not weights:
            weights = get_weights(pim1, pi)
        mutual_element = (1/weights)*analysis.py_rod_math.normalize(pim1) + (1/(1-weights))*analysis.py_rod_math.normalize(pi)
        return analysis.py_rod_math.normalize(mutual_element)

    def get_mutual_axes_inverse(mim1, mi, pim1=None, pi=None, weights=None):
        if not weights:
            weights = get_weights(pim1, pi)
        m_mutual = (mim1*(1.0/weights) + mi*(1.0/(1-weights)))/(analysis.py_rod_math.get_length(mi)+analysis.py_rod_math.get_length(mim1))
        return analysis.py_rod_math.normalize(m_mutual)
    
    try:
        analysis.p_i
    except AttributeError:
        analysis.p_i = analysis.get_p_i(analysis.rod.current_r)
    try:
        analysis.equil_p_ir
    except AttributeError:
        analysis.equil_p_i = analysis.get_p_i(analysis.rod.equil_r)
        
    #try:
    #    analysis.get_constant_EI()
    #except ValueError:
    #    import warnings
    #    warnings.warn("EI is not constant. If this is for an equipartition test, it won't work.")
        
    delta_omega = np.zeros( [len(analysis.p_i), len(analysis.p_i[0])-1, 2] )
    L_i = np.zeros( [len(analysis.p_i), len(analysis.p_i[0])-1] )

    
    for frame in range(len(analysis.p_i)):
        for element in range(len(delta_omega[0])):
            
            # NOTE TO SELF: ADD NORMALISATION!
            
            #B_j = np.matrix([[self.rod.B_matrix[frame][element][0], self.rod.B_matrix[0][0][1]], [self.rod.B_matrix[0][0][2], self.rod.B_matrix[0][0][3]]])
            #B= np.matrix([[analysis.rod.B_matrix[frame][element+1][0], analysis.rod.B_matrix[0][0][1]], [analysis.rod.B_matrix[0][0][2], analysis.rod.B_matrix[0][0][3]]])
            
            #mutual_l = get_mutual_frame_inverse(analysis.p_i[frame][element], analysis.p_i[frame][element+1])
            #equil_mutual_l = get_mutual_frame_inverse(analysis.equil_p_i[frame][element], analysis.equil_p_i[frame][element+1])
            
            pim1 = analysis.p_i[frame][element]
            pi = analysis.p_i[frame][element+1]
            mim1 = analysis.py_rod_math.normalize(analysis.rod.current_m[frame][element])
#            nim1 = np.cross(mim1, pim1/np.linalg.norm(pim1))
            mi = analysis.py_rod_math.normalize(analysis.rod.current_m[frame][element+1])
#            ni = np.cross(mi, pi/np.linalg.norm(pi))
        
            equil_pim1 = analysis.equil_p_i[frame][element]
            equil_pi = analysis.equil_p_i[frame][element+1]
            equil_mim1 = analysis.py_rod_math.normalize(analysis.rod.equil_m[frame][element])
#            equil_nim1 = np.cross(equil_mim1, equil_pim1/np.linalg.norm(equil_pim1))
            equil_mi = analysis.py_rod_math.normalize(analysis.rod.equil_m[frame][element+1])
#            equil_ni = np.cross(equil_mi, equil_pi/np.linalg.norm(equil_pi))
            
            if fast:
                energy, omega, omega_equil = FFEA_rod.rod_math.get_bend_energy(pim1, pi, equil_pim1, equil_pi, mim1, equil_mim1, mi, equil_mi, np.array([1,0,0,1]), np.array([1,0,0,1]) )
                delta_omega[frame][element] = omega - omega_equil
                L_i[frame][element] = (analysis.py_rod_math.get_length(equil_pi)+analysis.py_rod_math.get_length(equil_pim1))/2.0
                continue
            
            #get weighting between elements adjacent to bending node!
            weight = get_weights(pim1, pi)
            
            #get weighted element at that node
            mutual_l = get_mutual_frame_inverse(pim1, pi, weights=weight)

            # parallel transport our two existing material axes onto our new element
            mutual_mi = np.array(analysis.py_rod_math.parallel_transport(mi,  analysis.py_rod_math.normalize(pi), mutual_l))[0]
            mutual_mim1 = np.array(analysis.py_rod_math.parallel_transport(mim1,  analysis.py_rod_math.normalize(pim1), mutual_l))[0]
            
            # get the weighted average to create the mutual material axes
            mutual_m_rotated = get_mutual_axes_inverse(mutual_mim1, mutual_mi, weights=weight)
            
            omega, kb = analysis.py_rod_math.omega(pi, pim1, np.cross(mutual_l, mutual_m_rotated), mutual_m_rotated)
            
            #get weighting between elements adjacent to bending node!
            equil_weight = get_weights(equil_pim1, equil_pi)
            
            #get weighted element at that node
            equil_mutual_l = get_mutual_frame_inverse(equil_pim1, equil_pi, weights=equil_weight)

            # parallel transport our two existing material axes onto our new element
            equil_mutual_mi = np.array(analysis.py_rod_math.parallel_transport(equil_mi, analysis.py_rod_math.normalize(equil_pi), equil_mutual_l))[0]
            equil_mutual_mim1 = np.array(analysis.py_rod_math.parallel_transport(equil_mim1, analysis.py_rod_math.normalize(equil_pim1), equil_mutual_l))[0]
            
            # get the weighted average to create the mutual material axes
            equil_mutual_m_rotated = get_mutual_axes_inverse(equil_mutual_mim1, equil_mutual_mi, weights=equil_weight)
            
            equil_omega, equil_kb = analysis.py_rod_math.omega(equil_pi, equil_pim1, np.cross(equil_mutual_l, equil_mutual_m_rotated), equil_mutual_m_rotated)
            
            delta_omega[frame][element] = np.array(omega - equil_omega).transpose()[0]
            #delta_omega[frame][element] = omega - equil_omega
            
            L_i[frame][element] = (analysis.py_rod_math.get_length(equil_pi)+analysis.py_rod_math.get_length(equil_pim1))/2.0
            
    # return np.average(delta_omega, 0), np.average(L_i, 0)
    
    return delta_omega, L_i
            
            #inner = np.dot(np.transpose(delta_omega), np.dot(B, delta_omega))

            #energy = 0.5*(inner)*(1/((rod_math.get_length(equil_pi)+rod_math.get_length(equil_pim1))))
            
            
            
            #self.bending_energy[frame][element] = energy
            
def get_avg_and_mean_error(rod_energy, half_kBT, per_element = False):
    """
    For the rod energies, compute the average and standard error.
    Parameters:
        rod_energy: a 2-d array containing the energy for each node in each
        frame.
        half_hBT: boltzmann's constant * temperature * 0.5. A float.
        per_element: if False, will return the average energy. If True, the
        function will return an array of energies, one for each element.
    """
    if per_element:
        quantity=half_kBT/np.nanmean(rod_energy, axis=0)
        standard_error = (half_kBT/np.std(rod_energy, axis=0))/np.sqrt(len(rod_energy))
    else:
        quantity = half_kBT/np.nanmean(rod_energy)
        standard_error = (half_kBT/np.std(rod_energy))/np.sqrt(len(rod_energy))
    return quantity, standard_error

def get_B_avg(delta_omega, temperature, li, element_no, covariance=True, get_C=False):
    """
    \f$ B = k_B T l_i \cdot C^{-1} \f$
    Where
    \f$ C_{ij} = \langle \Delta \omega_i \Delta \omega_j \rangle =
    \begin{bmatrix} \langle \Delta \omega_1 \cdot \Delta \omega_1 \rangle, \langle \Delta \omega_1 \cdot \Delta \omega_2 \rangle \\ \langle \Delta \omega_2 \cdot \Delta \omega_1 \rangle, \langle \Delta \omega_2 \cdot \Delta \omega_2 \rangle \end{bmatrix} \f$
    This expression lets us recover the different elemenents of $B$, the bending
    stiffness matrix, from the values of $\Delta \omega$, the material curvature,
    and also the temperature $T$, Boltzmann's constant $k_B$, and $l_i$, the
    absolute length of the elements on either side of the node$.
    Parameters:
        delta_omega, a 2-d array containing the value of delta_omega for each
        frame and node.
        Temperature, the temperature in K.
        element_no: the index of the element to get the B matrix for.
    Returns:
        B, a numpy array representing the B matrix.
    """
    #delta_omega = omega_tilde-omega
    
    delta_omega_element = delta_omega[:,element_no]
    li_element = li[:,element_no]
    if covariance:
        C11 = np.average(delta_omega_element[:,0]*delta_omega_element[:,0]) - (np.average(delta_omega_element[:,0])*np.average(delta_omega_element[:,0]))
        C12_21 = np.average(delta_omega_element[:,0]*delta_omega_element[:,1]) - (np.average(delta_omega_element[:,0])*np.average(delta_omega_element[:,1]))
        C22 = np.average(delta_omega_element[:,1]*delta_omega_element[:,1]) - (np.average(delta_omega_element[:,1])*np.average(delta_omega_element[:,1]))
    else:
        C11 = np.average(delta_omega_element[:,0]*delta_omega_element[:,0])
        C12_21 = np.average(delta_omega_element[:,0]*delta_omega_element[:,1])
        C22 = np.average(delta_omega_element[:,1]*delta_omega_element[:,1])
    C = np.matrix( [[C11, C12_21],[C12_21, C22]])
    try:
        B = (1.38064852e-23*temperature* (np.average(li_element))) * np.linalg.inv(C)
    except np.linalg.LinAlgError:
        print(C)
        raise np.linalg.LinAlgError()
    if get_C:
        return B, C
    return B

    # C from initial recovery is too small  (therefore B_old is over-estimate)
    # C for second run should be smaller, so C - C_old should be positive

def iterative_improve_B(original_rod, new_rod, write_out_file="", extra_data = True, alpha=1.0):
    """
    An iterative method to recover the value of B from two successive rod
    trajectores. The first trajectory is any trajectory. The second is a
    trajectory based on the values of B recovered from the first trajectory.
    The value of B is recovered using the following equation:
    \f[\bm{B}_{new} = k_B T \bm{L}_i \left( \bm{C}_i - \bm{C^\prime}_i + k_B T \bm{L_i} \bm{B}_i^{-1} \right)^{-1}\f]
    Where \f$\bm{B}_{new}\f$ is the corrected value of \f$\bm{B}$, $k_B\f$ is
    Boltzmann's constant, \f$T\f$ is the temperature, \f$\bm{L}_i\f$ is
    \f$L_i\f$, defined in the method FFEA_rod.get_L_i, \f$\bm{C}_i\f$ is the
    \f$\bm{C}\f$ matrix for the original trajectory \f$\Gamma\f$,
    \f$\bm{C^\prime}_i\f$ is the \f$\bm{C}\f$ matrix for the new trajectory
    \f$\Gamma^\prime\f$, and \f$\bm{B}_i\f$ is the bending stiffness matrix
    calculated for the original trajectory.
    
    Params:
        original_rod - FFEA_rod object, the rod trajectory of the rod Gamma.
        new_rod - FFEA_rod object, rod trajectory of the rod Gamma'.
        write_out_file - string, if this is non-empty, a new rod file will be
        written to that path with the new values of B
        extra_data - bool, only for debugging. See return values.
        
    Returns:
        if extra_data is false: B, as a numpy array
        if extra_data: B, old_B, old_C, new_B, new_C
    
    """
    temp = 300 #kelvin
    analytical_kbT = temp*1.38064852 * 10**-23
    half_kbT = 0.5*analytical_kbT
    
    old_analysis = FFEA_rod.anal_rod(original_rod)
    new_analysis = FFEA_rod.anal_rod(new_rod)
    
    old_delta_omega, old_L_i = get_delta_omega(old_analysis, fast=True)
    old_B = []
    old_C = []
    for element_no in range(len(old_delta_omega[0])):
        curr_old_B, curr_old_C = get_B_avg(old_delta_omega, temp, old_L_i, element_no, get_C = True, covariance=True)
        old_B.append(curr_old_B)
        old_C.append(curr_old_C)
        
    new_delta_omega, new_L_i = get_delta_omega(new_analysis, fast=True)
    new_B = []
    new_C = []
    for element_no in range(len(new_delta_omega[0])):
        curr_new_B, curr_new_C = get_B_avg(new_delta_omega, temp, new_L_i, element_no, get_C = True, covariance=True)
        new_B.append(curr_new_B)
        new_C.append(curr_new_C)
    
    corrected_B = []
    #corrected_B2 = []
    for i in range(len(new_B)):
        corrected_B.append( analytical_kbT * np.average(new_L_i, 0)[i] * np.linalg.inv( alpha*(old_C[i] - new_C[i]) + (analytical_kbT * np.average(new_L_i, 0)[i] * np.linalg.inv( old_B[i] )) ) )
        #corrected_B2.append(  np.linalg.inv( (1/(analytical_kbT*np.average(new_L_i, 0)[i])*( old_C[i] - new_C[i] )) + np.linalg.inv(old_B[i])  ) )

    if write_out_file != "":
        fixed_rod = copy.copy(new_rod)
        fixed_rod.num_frames = 1
        fixed_rod.B_matrix[0][1:-1] = np.reshape(np.array(corrected_B).flatten(), [len(np.array(corrected_B).flatten())/4, 4])
        fixed_rod.write_rod(write_out_file)

    if extra_data:
        return np.reshape(np.array(corrected_B).flatten(), [len(np.array(corrected_B).flatten())/4, 4]), old_B, old_C, new_B, new_C, analytical_kbT * np.average(new_L_i, 0)

    return np.reshape(np.array(corrected_B).flatten(), [len(np.array(corrected_B).flatten())/4, 4])

def recover_B(rodtraj_path, new_rod_path=None, fast=True):
    """
    Recover only B from a given rod trajectory. Uses the method defined in
    get_B_avg.
    Params:
        - rodtraj_path - string, path to a rod trajectory file
        - new_rod_path -  string. If this is set, a new .rod file will be
        created with the B values recovered from the original trajectory.
        - fast - bool, if true, will try and use the fast math functions defined
        in rod_math_core. If rod_math_core failed to compile, this won't work.
    Returns:
        B, as a 2d numpy array. This would correspond to rod.B_matrix[0].
    """
    temp = 300
    print("Loading...")
    rod = FFEA_rod.FFEA_rod(rodtraj_path)
    analysis = FFEA_rod.anal_rod(rod)
    print("Recovering B...")
    delta_omega, L_i = get_delta_omega(analysis, fast=fast)
    B = []
    for element_no in range(len(delta_omega[0])):
        B.append(get_B_avg(delta_omega, temp, L_i, element_no))
    if new_rod_path:
        rod.num_frames = 1
        rod.B_matrix[0] = np.zeros( np.shape(rod.B_matrix[0]) )
        rod.B_matrix[0][1:-1] = np.reshape(np.array(B).flatten(), [len(np.array(B).flatten())/4, 4])
        rod.B_matrix[0][0] = rod.B_matrix[0][1]
        rod.B_matrix[0][-1] = rod.B_matrix[0][-2]
        rod.write_rod(new_rod_path)
    return np.array(B)

def complete_parallel_transport(analysis, from_index, to_index, equil=True, force_normalize=True):

    if force_normalize:
        n = FFEA_rod.rod_math.normalize
    else:
        n = lambda n: n
    
    try:
        analysis.p_i
    except AttributeError:
        analysis.p_i = analysis.get_p_i(analysis.rod.current_r)
        analysis.equil_p_i = analysis.get_p_i(analysis.rod.equil_r)
    
    if equil:
        m = analysis.rod.current_m
        p = analysis.equil_p_i
    else:
        m = analysis.rod.equil_m
        p = analysis.p_i
    
    curr_mataxis = n(m[0][from_index])
    for i in range(from_index,to_index-1):
        curr_mataxis = np.asarray(FFEA_rod.rod_math.parallel_transport(curr_mataxis, n(p[0][i]), n(p[0][i+1])))[0]
    return curr_mataxis

def get_material_axis_angle(analysis, node_index, frame=0, equil=True):
    
    n = FFEA_rod.rod_math.normalize
    
    try:
        analysis.p_i
    except AttributeError:
        analysis.p_i = analysis.get_p_i(analysis.rod.current_r)
        analysis.equil_p_i = analysis.get_p_i(analysis.rod.equil_r)
    
    m1 = analysis.rod.equil_m[0][node_index]
    m2 = analysis.rod.equil_m[0][node_index+1]
    p1 = n(analysis.equil_p_i[0][node_index])
    p2 = n(analysis.equil_p_i[0][node_index+1])
    
    m1_prime = FFEA_rod.rod_math.parallel_transport(m1, p1, p2)
    angle = FFEA_rod.rod_math.get_signed_angle(m1_prime, m2, p2)
    return angle

def unroll(analysis, perpendicularize=True, do_complete_parallel_transport=False):
    print("Unrolling material axes...")
    try:
        analysis.p_i
    except AttributeError:
        print("Generating p_i...")
        analysis.p_i = analysis.get_p_i(analysis.rod.current_r)
        analysis.equil_p_i = analysis.get_p_i(analysis.rod.equil_r)
        analysis.p_i_equil = analysis.equil_p_i
        
    rod = analysis.rod
    
    print("Normalizing...")
    for frame in range(len(rod.current_m)):
        for elem in range(len(rod.current_m[frame])-1):
            rod.current_m[frame][elem] = FFEA_rod.rod_math.normalize(rod.current_m[frame][elem])
            rod.equil_m[frame][elem] = FFEA_rod.rod_math.normalize(rod.equil_m[frame][elem])
    
    print("Perpendicularizing...")
    if perpendicularize:
        for frame in range(len(rod.current_m)):
            for elem in range(len(rod.current_m[frame])-1):
                rod.current_m[frame][elem] = FFEA_rod.rod_math.perpendicularize(rod.current_m[frame][elem], analysis.p_i[frame][elem])
                rod.equil_m[frame][elem] = FFEA_rod.rod_math.perpendicularize(rod.equil_m[frame][elem], analysis.p_i_equil[frame][elem])
    
    # some init stuff
    
    rotation_table = np.zeros(len(rod.current_m[0]))
    mat1 = FFEA_rod.rod_math.normalize(rod.equil_m[0][0])
    elem1 = FFEA_rod.rod_math.normalize(analysis.p_i_equil[0][0])
    
    print("Creating rotation table...")
    for mataxis_index in range(len(rod.equil_m[0])-1):
        if do_complete_parallel_transport:
            transported_mataxis = complete_parallel_transport(analysis, 0, mataxis_index, equil=True, force_normalize=True)
            curr_mataxis = FFEA_rod.rod_math.normalize(rod.equil_m[0][mataxis_index])
            curr_elem = FFEA_rod.rod_math.normalize(analysis.p_i_equil[0][mataxis_index])
            rotation_table[mataxis_index] = FFEA_rod.rod_math.get_signed_angle(transported_mataxis, curr_mataxis, curr_elem)
            
        else:
            transported_mataxis = np.squeeze(np.asarray(FFEA_rod.rod_math.parallel_transport(FFEA_rod.rod_math.normalize(rod.equil_m[0][mataxis_index]), FFEA_rod.rod_math.normalize(analysis.p_i_equil[0][mataxis_index]), elem1)))
            # note: transported mataxis is np array type
            rotation_table[mataxis_index] = FFEA_rod.rod_math.get_signed_angle(mat1, transported_mataxis, elem1)

    print("Applying rotation table...")
    for frame in range(len(rod.current_m)):
        if frame%10000 == 0:
            print("Frame "+str(frame)+" of "+str(rod.num_frames))
        for elem in range(len(rod.current_m[frame])-1):
            k = FFEA_rod.rod_math.normalize(analysis.p_i[frame][elem])
            rod.current_m[frame][elem] = FFEA_rod.rod_math.rodrigues(FFEA_rod.rod_math.normalize(rod.current_m[frame][elem]), k, rotation_table[elem])
            k_equil = FFEA_rod.rod_math.normalize(analysis.p_i_equil[frame][elem])
            rod.equil_m[frame][elem] = FFEA_rod.rod_math.rodrigues(FFEA_rod.rod_math.normalize(rod.equil_m[frame][elem]), k_equil, rotation_table[elem])

    analysis.rod = rod
    return analysis

#prmtop_file = "/home/rob/pCloudDrive/pCloud Sync/scratchpad/samsim/SMC_extract/short11799SMC_open_no_atpmd4.top"
#mdcrd_traj_file = "/home/rob/pCloudDrive/pCloud Sync/scratchpad/samsim/SMC_extract/short11799SMC_open_no_atpmd4.x"
#inpcrd_file = "/home/rob/pCloudDrive/pCloud Sync/scratchpad/samsim/SMC_extract/short11799SMC_open_no_atpmd4.rst7"
#get_B = True
#target_length=15
#cluster_size = 10
#simplify_only = False
#rod=None
#get_inhomogeneous_beta=True
#get_inhomogeneous_kappa=True
#unroll_rod = True
#rod_out = "/home/rob/pCloudDrive/pCloud Sync/scratchpad/samsim/SMC_extract/SMC.rodtraj" # CHANGE ME

# out = main(prmtop_file, mdcrd_traj_file, inpcrd_file, get_B = True, target_length=15, cluster_size = 10, simplify_only = False, rod=None, radius = 5e-9, rod_out=rod_out, unroll_rod=True, get_inhomogenous_beta=True, get_inhomogenous_kappa=True)

def main(prmtop_file, mdcrd_traj_file, inpcrd_file, chain1_end, trj_format="TRJ", inp_format="INPCRD", get_B = True, target_length=15, cluster_size = 10, simplify_only = False, rod=None, radius = 5e-9, rod_out=None, unroll_rod=True, get_inhomogenous_beta=True, get_inhomogenous_kappa=True):
    """
    For an atomistic trajectory, this will create an equivalent rod trajectory,
    and use that trajectory to compute the material parameters of the rod.
    Cool huh?
    Parameters:
        prmtop_file, mdcrd_tarj_file, inpcrd_file: paths to AMBER output files
        (strings)
        chain1_end: the index of the atom at which the first chain (the first
        coil) in the coiled-coil ends (int)
        trj_format and inp_format: formats of the trajectory and equilibrium
        structure files (strings). Uses MDAnalysis formats.
        get_B: whether to also compute the anisotropic B matrices (bool)
        target_length: the number of elements in the rod to be created (int)
        simplify_only: if True, don't parameterise rod, just make it (bool)
        rod: if you already have a rod trajectory, supply a rod object here.
    Returns:
        analysis, an instance of FFEA_rod.anal_rod
        delta_omega, L_i, B - the values used to compute the anisotropic B
        matrix. Arrays.
    It also prints out the isotropic results, for reference.
    """
    if not rod:
        
        print("Loading files...")
        topology = prmtop_file
        trajectory = mdcrd_traj_file
        
        u_initial = MDAnalysis.Universe(topology, inpcrd_file, format=inp_format)
        initial_backbone = u_initial.select_atoms('protein and backbone').positions/1e10
        #initial_chain_len = (np.shape(initial_backbone)[0])/2
        #initial_chain2 = initial_backbone[:initial_chain_len]
        #initial_chain1 = initial_backbone[initial_chain_len:]
        #initial_r = (initial_chain1+initial_chain2)/2
        #initial_m = initial_chain2 - initial_chain1
        
        u = MDAnalysis.Universe(topology, trajectory, format=trj_format)
        backbone = u.select_atoms('protein and backbone')
        backbone_traj = np.zeros([len(u.trajectory), len(backbone), 3])
        
        print("Retrieving backbone trajectory (warning: slow as hell, blame MDAnalysis)...")
        for frame_no in range(len(u.trajectory)):
            backbone_traj[frame_no] = backbone.positions
            try:
                u.trajectory.next()
            except ValueError:
                break
            except StopIteration:
                break
        
        #chain1_end = 1372
#        chain2_hinge = backbone_traj[:][1997:2023]
#        chain1_hinge_len = 78
#        chain2_hinge_len = len(chain2_hinge)
        #pivot_index = 624
        
        print("Constructing rod...")
        
        chain_len = (np.shape(backbone_traj)[1])/2
        
        backbone_traj /= 1e10 # convert from A to m

        #1) work out how many elements we need, create a rod
        #2) set the indices of averaging clusters
        #3) get the nearest (opposite) element lookup table
        
        def set_cluster_indices(target_length, cluster_size, chain1_end):
            """
            To construct the rod trajectory, only small clusters of atoms are
            averaged. This function sets the indexes of the nodes which will
            be used in each cluster.
            Params:
                target_length: target number of nodes, an int
                cluster_size: width (in atoms) of each cluster
                chain1_end: index of the atom at the end of the first chain
            Returns:
                cluster_indices, a 2-d array of integers where the first axis
                is the rod node id and the second is the atomistic atom id.
            """
            # note: these indices are defined for the FIRST chain only.
            cluster_indices = np.zeros([target_length,cluster_size], dtype=int)
            starting_cluster_index = int(cluster_size/2)
            end_cluster_index = chain1_end - int(cluster_size/2)
            cluster_center_indices = np.linspace(starting_cluster_index, end_cluster_index, target_length, dtype=int)

            for cluster_index in range(len(cluster_indices)):
                i = 0
                for node_index in range(len(cluster_indices[cluster_index])):
                    cluster_indices[cluster_index][node_index] = cluster_center_indices[cluster_index]+i
                    if i == 0:
                        i+=1
                    elif i>0:
                        i *= -1
                    elif i<0:
                        i *= -1
                        i += 1
            return cluster_indices
        
        def get_nearest_node_table(initial_backbone, chain1_end):
            """
            Get a lookup table of the index of the nearest node in the opposite
            chain for each node.
            Params:
                initial_backbone - the initial (equilibrium) state of the
                atomistic structure.
                chain1_end - integer, index of the atom at the end of the first
                chain
            Returns:
                nearest_node_lookup_table, a 1-d array of integers where the
                index is the node index and the value is the index of the
                nearest node in the opposite chain.
            """
            chain1 = initial_backbone[:chain1_end]
            chain2 = initial_backbone[chain1_end:]
            chain1_nearest_node_table = np.zeros(len(chain1), dtype=int)
            #chain2_nearest_node_table = np.zeros(len(chain2), dtype=int)
            for node_index in range(len(chain1)):
                distances_vec = chain1[node_index] - chain2
                distances_scalar = np.linalg.norm(distances_vec, axis=1)
                min_index = np.argmin(distances_scalar)
                chain1_nearest_node_table[node_index] = min_index
                #chain2_nearest_node_table[min_index] = node_index
            #return np.concatenate([chain1_nearest_node_table, chain2_nearest_node_table])
            return chain1_nearest_node_table+chain1_end
            
        def get_node_traj(cluster_indices, node_index, backbone_traj, nearest_node_table, chain1_end):
            """
            Given the cluster indices and nearest nodes, (see above), get the
            cluster-averaged position of a particular rod node, for the whole
            trajectory. The position is the average of 10 atoms in both chains.
            The material axis is the average vector that goes between those
            atoms (normalized).
            Params:
                cluster_indices: ndarray generated by get_cluster_indices
                node_index: index of the rod node to find a trajectory for
                backbone_traj: all-atom trajectory from mdanalysis as an ndarray
                nearest_node_table: generated from get_nearest_node_table
            Returns:
                node position as a http://www.np.array/3-d array where the first index is the frame,
                second axis is the node index, and third axis is the dimension.
                material axis vector in the same structure.
            """
            cluster_size = len(cluster_indices[0])
            chain1_nodes = np.zeros([cluster_size, len(backbone_traj), 3])
            chain2_nodes = np.zeros([cluster_size, len(backbone_traj), 3])
            i=0
            for curr_node_index in cluster_indices[node_index]:
                chain1_nodes[i] = backbone_traj[:,curr_node_index]
                chain2_nodes[i] = backbone_traj[:,nearest_node_table[curr_node_index-1]]
                i+=1
            chain1_avg = np.average(chain1_nodes, axis=0)
            chain2_avg = np.average(chain2_nodes, axis=0)
            
            #setup mataxis_cluster
            mataxis_chain1_nodes = np.zeros([cluster_size, len(backbone_traj), 3])
            mataxis_chain2_nodes = np.zeros([cluster_size, len(backbone_traj), 3])
            mataxis_offset = int(np.average(cluster_indices[1] - cluster_indices[0])/2.0)
            
            j=0
            for curr_node_index in cluster_indices[node_index]:
                if curr_node_index+mataxis_offset > chain1_end:
                    break # we're at the end of the rod and this is error handling for control flow becauseit's 6AM and I'm sad
                mataxis_chain1_nodes[j] = backbone_traj[:,curr_node_index+mataxis_offset]
                mataxis_chain2_nodes[j] = backbone_traj[:,nearest_node_table[curr_node_index+mataxis_offset-1]]
                j+=1
            
            mataxis = mataxis_chain1_nodes - mataxis_chain2_nodes
            m_average = np.average(mataxis, axis=0)
            absolutes = np.linalg.norm(m_average, axis=1)
            absolutes = np.array([absolutes, absolutes, absolutes])
            absolutes = np.swapaxes(absolutes, 0, 1)
            m_average = m_average/absolutes
            
            return np.average([chain1_avg, chain2_avg], axis=0), m_average
        
        print("Setting cluster indices...")
        cluster_indices = set_cluster_indices(target_length, cluster_size, chain1_end)
        print("Calculating nearest node table...")
        nearest_node_table = get_nearest_node_table(initial_backbone, chain1_end)
        print("Calculating node positions and material axes...")
        rod_r = np.zeros([len(backbone_traj), target_length, 3])
        rod_m = np.zeros([len(backbone_traj), target_length, 3])
        for node_index in range(target_length):
            rod_r[:,node_index], rod_m[:,node_index] = get_node_traj(cluster_indices, node_index, backbone_traj, nearest_node_table, chain1_end)

        equil_r = np.zeros([len(backbone_traj), target_length, 3])
        equil_m = np.zeros([len(backbone_traj), target_length, 3])
        for node_index in range(target_length):
            equil_r[:,node_index], equil_m[:,node_index] = get_node_traj(cluster_indices, node_index, np.array([initial_backbone]), nearest_node_table, chain1_end)
   
        rod_m[np.isnan(rod_m)] = 0
        equil_m[np.isnan(equil_m)] = 0 # old trick
    
        print("Initializing FFEA rod object...")
        rod = FFEA_rod.FFEA_rod(num_elements=target_length)
        rod.current_r = rod_r
        rod.current_m = rod_m
        
        rod.equil_r = equil_r
        rod.equil_m = equil_m
        
        rod.material_params = np.ones([len(backbone_traj), target_length, 3])
        rod.B_matrix = np.ones([len(backbone_traj), target_length, 4])
        
        rod.perturbed_x_energy_positive = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.perturbed_y_energy_positive = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.perturbed_z_energy_positive = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.perturbed_x_energy_negative = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.perturbed_y_energy_negative = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.perturbed_z_energy_negative = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.twisted_energy_positive = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
        rod.twisted_energy_negative = np.zeros([len(rod.current_m), len(rod.current_m[0]), 3])
    
    print("Setting parameters...")
    
    
    FFEA_rod.rod_creator.set_params(rod, 1, 1, radius, 1)
    rod.num_frames = len(rod.current_r)
    
    if simplify_only:
        return rod
    
    analysis = FFEA_rod.anal_rod(rod)
    
    if unroll_rod:
        print("Unrolling...")
        analysis = unroll(analysis)
        
    try:
        analysis.p_i
        analysis.equil_p_i
    except AttributeError:
        rod.p_i = rod.get_p_i(rod.current_r)
        rod.equil_p_i = rod.get_p_i(rod.equil_r)
        analysis.p_i = analysis.rod.p_i
        analysis.p_i_equil = analysis.rod.equil_p_i
        analysis.equil_p_i = analysis.rod.equil_p_i
    
    temp = 300 #kelvin
    analytical_kbT = temp*1.38064852 * 10**-23
    
    half_kbT = 0.5*analytical_kbT
    
    print("Computing average constants...")
    analysis.get_stretch_energy()
    half_stretch_squared_avg = np.average(analysis.stretch_energy)
    kappa_stretch = half_kbT/half_stretch_squared_avg
    
    analysis.get_bending_response_mutual()
    half_bend_squared_avg = np.nanmean(analysis.bending_energy)
    B_isotropic = half_kbT/half_bend_squared_avg
    
    analysis.get_twist_amount()
    half_twist_squared_avg = np.average(analysis.twist_energy)
    beta = half_kbT/half_twist_squared_avg
    
    avg_stretch, stretch_error = get_avg_and_mean_error(analysis.stretch_energy, half_kbT)
    avg_bend, bend_error = get_avg_and_mean_error(analysis.bending_energy, half_kbT)
    avg_twist, twist_error = get_avg_and_mean_error(analysis.twist_energy, half_kbT)
    
    print("Done.")
    print("Kappa (stretch constant)    = "+str(kappa_stretch)+"+/-"+str(stretch_error))
    print("B (isotropic bend constant) = "+str(B_isotropic)+"+/-"+str(bend_error))
    print("beta (twist constant)       = "+str(beta)+"+/-"+str(twist_error))
    
    if not get_inhomogenous_beta and not get_inhomogenous_kappa and not get_B:
        return analysis
    
    if get_inhomogenous_beta or get_inhomogenous_kappa:
        print("Computing inhomogeneous constants...")
    
    if get_inhomogenous_beta:
        inhomogeneous_beta, beta_err = get_inhomogeneous_param(analysis, 1, kbT=analytical_kbT)
        inhomogeneous_beta = np.concatenate([[inhomogeneous_beta[0]], inhomogeneous_beta, [inhomogeneous_beta[-1]] ])
        analysis.rod.material_params[:,:,1] = inhomogeneous_beta
        print("beta (inhomogeneous)"+str(inhomogeneous_beta))
        
    if get_inhomogenous_kappa:
        inhomogeneous_kappa, kappa_err = get_inhomogeneous_param(analysis, 0, kbT=analytical_kbT)
        inhomogeneous_kappa = np.concatenate([inhomogeneous_kappa, [inhomogeneous_kappa[-1]] ])
        analysis.rod.material_params[:,:,0] = inhomogeneous_kappa
        print("kappa (inhomogeneous)"+str(inhomogeneous_kappa))
    
    if get_B:
        print("Computing B...")
        delta_omega, L_i = get_delta_omega(analysis, fast=False)
        B = []
        for element_no in range(len(delta_omega[0])):
            B.append(get_B_avg(delta_omega, temp, L_i, element_no))

        B_flat = np.ndarray.flatten(np.asarray(B))
        B_arr = np.zeros([np.shape(delta_omega)[1]+2, 4 ])
        B_arr[1:-1] = np.reshape(B_flat, [np.shape(delta_omega)[1], 4])
        B_arr[0] = B_arr[1]
        B_arr[-1] = B_arr[-2]
        analysis.rod.B_matrix[:] = B_arr
        print("Saving rod...")
        print("B (inhomogeneous, anisotropic): "+str(B_arr))
        print("Writing rod...")

    if rod_out:
        analysis.rod.write_rod(rod_out)

    return analysis, delta_omega, L_i, np.array(B), inhomogeneous_beta, inhomogeneous_kappa
    
#if __name__ == '__main__' and '__file__' in globals():
#    args = parser.parse_args()
#    # Get args and build objects
#    out = main(args.prmtop_file, args.mdcrd_file, args.inpcrd_file, args.chain1_end, args.trj_format, args.inp_format, get_B = args.get_B, target_length=args.target_length, cluster_size = args.cluster_size, simplify_only = False, rod=None, radius = args.radius, rod_out=args.rod_out, unroll_rod=args.unroll_rod, get_inhomogenous_beta=args.get_inhomogeneous_beta, get_inhomogenous_kappa=args.get_inhomogeneous_kappa)
#    print out
