# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 21:10:35 2021

@author: z003p2nh
"""
from scipy.interpolate import CubicHermiteSpline as hspline
import numpy as np

# constants
z_max = 0.389
z_op = 0.25    # usual operating hight
z_s = -0.02    # step hight
T_c = 1.5      # cycle time

<<<<<<< HEAD
#phase = np.array([0.75, 0.5, 0.25, 0.])
phase = np.array([0.0, 0.5, 0.0, 0.5])
=======
phase = np.array([0.75, 0.5, 0.25, 0.])
>>>>>>> quaroDance/main
support_ratio = 0.8
stride = 0.05
swing_spline = hspline([0, 0.5, 1],
                       [0, z_s, 0],
                       [0,   0, 0])
'''
# helping lists
tsteps = [(n/8.)*T_c for n in range(8+1)]

# spline generation
get_z_trot = hspline(tsteps,
                     [[0,z_s,0,0,0,0,0,0,0],
                      [0,0,0,z_s,0,0,0,0,0],
                      [0,0,0,0,0,z_s,0,0,0],
                      [0,0,0,0,0,0,0,z_s,0]], 
                     [[0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0],
                      [0,0,0,0,0,0,0,0,0]],
                     axis=1, extrapolate='periodic')
'''



def get_x_trot(t):
    '''return x position of each leg in local coordinates'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    supporting = t_phi <= support_ratio
    x = np.zeros(4)
    x_stance = stride/2 - stride*(t_phi/support_ratio)
    x_swing = -stride/2 + stride*((t_phi - support_ratio)/(1. - support_ratio))
    x = x_swing
    x[supporting] = x_stance[supporting]
    return x.T

def get_y_trot(t):
    '''return y position of each leg in local coordinates'''
    return np.zeros((4, len(t)))

def get_z_trot(t):
    '''return z position of each leg in local coordinates'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    supporting = t_phi <= support_ratio
    z = np.zeros(t_mat.T.shape)
    z_swing = swing_spline((t_phi - support_ratio)/(1. - support_ratio))
    z[~supporting] = z_swing[~supporting]
    return z.T

def get_legstate_trot(t):
    '''return legstate as bool'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    return t_phi <= support_ratio


'''
get_z_jump = hspline([0,       2,    2.05,    3,    4], 
                [z_op, 0.25,   z_max, z_op, z_op], 
                [0,       0,       0,    0,    0])
'''
