# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 21:10:35 2021

@author: z003p2nh
"""
from scipy.interpolate import CubicHermiteSpline as hspline
import numpy as np
from . import zero_moment_point as zmp

# constants
z_max = 0.389
z_op = 0.35    # usual operating hight
z_s = -0.02    # step hight
T_c = 1      # cycle time
N_spline = 30  # number of reference points for ZMP spline


#phase = np.array([0.75, 0.5, 0.25, 0.])
phase = np.array([0., 0.5, 0.5, 0.])
#phase = np.array([0.0, 0.5, 0.0, 0.5])
support_ratio = 0.8
stride = 0.15
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
    '''return x leg position vector in local coordinates'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    supporting = t_phi <= support_ratio
    x = np.zeros(4)
    x_stance = stride/2 - stride*(t_phi/support_ratio)
    x_swing = -stride/2 + stride*((t_phi - support_ratio)/(1. - support_ratio))
    x = x_swing
    x[supporting] = x_stance[supporting]
    return x

def get_y_trot(t):
    '''return y leg position vector in local coordinates'''
    if (type(t) == int) or (type(t) == float):
        y =  np.zeros(4)
    else:
        y = np.zeros((len(t), 4))
    return y

def get_z_trot(t):
    '''return z leg position vector in local coordinates'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    supporting = t_phi <= support_ratio
    z = np.zeros(t_mat.T.shape)
    z_swing = swing_spline((t_phi - support_ratio)/(1. - support_ratio))
    z[~supporting] = z_swing[~supporting]
    return z

def get_legstate_trot(t):
    '''return legstate as bool'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    return t_phi <= support_ratio

def get_optimized_trot(t):
    '''returns 3x4 leg coordinates considering optimized leg position'''
    
    # calculate ZMP
    t_ref = np.linspace(0, t[-1], N_spline)
    x_ref = get_x_trot(t_ref).T
    y_ref = get_y_trot(t_ref).T
    leg_state = get_legstate_trot(t_ref)
    p_ref = np.zeros((2, N_spline))
    for n in range(N_spline):
        p_ref[0, n] = np.sum(x_ref[leg_state[n], n])/np.sum(leg_state[n])
        p_ref[1, n] = np.sum(y_ref[leg_state[n], n])/np.sum(leg_state[n])
        
    # calculate ZMP spline
    spline_zmp = zmp.optimize_ZMP_spline(p_ref, t[-1], 
                                         alpha=1e4, beta=1e-4, gamma=1e-6)
    
    # calculate desired leg points
    p_zmp = spline_zmp(t)
    x = get_x_trot(t) - np.tile(p_zmp[0, np.newaxis].T, 4)
    y = get_y_trot(t) - np.tile(p_zmp[1, np.newaxis].T, 4)
    z = get_z_trot(t) + z_op
    
    ### JUST A TEST; VECTORIZE THIS!
    N = len(t)
    pos = np.zeros((N, 3, 4))
    for n in range(N):
        pos[n,0,:] = x[n] 
        pos[n,1,:] = y[n] 
        pos[n,2,:] = z[n] 
    return pos


'''
get_z_jump = hspline([0,       2,    2.05,    3,    4], 
                [z_op, 0.25,   z_max, z_op, z_op], 
                [0,       0,       0,    0,    0])
'''
