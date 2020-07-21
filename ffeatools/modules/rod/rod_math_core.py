#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 12:03:37 2019

@author: rob
"""

inlinec = """
int main(){
    printf('Hello');
    return 0;
}
"""

import os
import subprocess
out = subprocess.check_output("gcc -shared -fPIC -O3 -o "+os.path.dirname(os.path.realpath(__file__))+os.path.sep+"librodmath.so "+os.path.dirname(os.path.realpath(__file__))+os.path.sep+"rod_math_core.c", shell=True)
dir_path = os.path.dirname(os.path.realpath(__file__))
#os.system("echo '"+inlinec+"' | gcc -x -shared -o librodmath.so -O3")
#os.system("echo -e '"+inlinec.replace("\n", "")+"' | gcc -xc -O3 -shared -o inline.so ")

import numpy as np
import ctypes
ndpointer = np.ctypeslib.ndpointer
rod_math_c = ctypes.cdll.LoadLibrary(os.path.join(dir_path, "librodmath.so"))

def arr(*args, **kwargs):
    kwargs.setdefault("dtype", np.float32)
    return np.array(*args, **kwargs)

cint = ctypes.c_int
cfloat = ctypes.c_float

rod_math_c.get_stretch_energy.restype = cfloat
rod_math_c.get_stretch_energy.argtypes = [cfloat, ndpointer(cfloat, flags="C_CONTIGUOUS"), ndpointer(cfloat, flags="C_CONTIGUOUS"),]

rod_math_c.get_bend_energy_mutual_parallel_transport.restype = cfloat
rod_math_c.get_bend_energy_mutual_parallel_transport.argtypes = [
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ]

rod_math_c.get_twist_energy.restype = cfloat
rod_math_c.get_twist_energy.argtypes = [
    cfloat,
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ndpointer(cfloat, flags="C_CONTIGUOUS"),
    ]

get_twist_energy = rod_math_c.get_twist_energy
get_bend_energy = rod_math_c.get_bend_energy_mutual_parallel_transport
get_stretch_energy = rod_math_c.get_stretch_energy

def get_twist_32(m_im1, m_im1_equil, m_i, m_i_equil, p_im1, p_im1_equil, p_i, p_i_equil, beta):
    return get_twist_energy(beta, m_i.astype("float32"), m_im1.astype("float32"), m_i_equil.astype("float32"),m_im1_equil.astype("float32"),p_im1.astype("float32"),p_i.astype("float32"),p_im1_equil.astype("float32"),p_i_equil.astype("float32"))

def get_stretch_32(k, p_i, p_i_equil):
    return get_stretch_energy( k, p_i.astype("float32"), p_i_equil.astype("float32") )

def get_bend_32(
        p_im1,
        p_i,
        p_im1_equil,
        p_i_equil,
        m_im1,
        m_im1_equil,
        m_i,
        m_i_equil,
        B_i_equil,
        B_im1_equil):
    omega = np.empty(2, dtype="float32")
    omega_equil = np.empty(2, dtype="float32")
    return get_bend_energy(p_im1.astype("float32"), p_i.astype("float32"), p_im1_equil.astype("float32"), p_i_equil.astype("float32"), m_im1.astype("float32"), m_im1_equil.astype("float32"), m_i.astype("float32"), m_i_equil.astype("float32"), B_i_equil.astype("float32"), B_im1_equil.astype("float32"), omega, omega_equil), omega, omega_equil