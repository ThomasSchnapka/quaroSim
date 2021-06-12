"""
Custom leg trajectories can be set up here. There are two main way, to create
a leg trajectory: using CubicHermiteSpline or an gait generation algorithm.

After creating the spline, it makes sense to use the zero_moment_point module
to create leg movements that follow the desired leg trajectory but still keep
the robot stable

"""
from scipy.interpolate import CubicHermiteSpline as hspline
import numpy as np
try:
    from . import zero_moment_point as zmp
except:
    pass

# constants
z_max = 0.389
z_op = 0.35    # usual operating hight
z_s = -0.05    # step hight
T_c = 1      # cycle time
N_ref = 30  # number of reference points for ZMP spline


#phase = np.array([0.75, 0.5, 0.25, 0.])    # walk
#phase = np.array([0., 0.5, 0.5, 0.])       # trot
phase = np.array([0.0, 0.5, 0.0, 0.5])     # gallop
#phase = np.array([0.0, 0.0, 0.5, 0.5])      # pronk
support_ratio = 0.8
stride = 0.1
swing_spline = hspline([0, 0.5, 1],
                       [0, z_s, 0],
                       [0,   0, 0])
leg_pos = np.array([[0.2, 0.2, -0.2, -0.2],
                    [0.15, -0.15, 0.15, -0.15],
                    [  0,    0,   0,    0]])
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



def get_x_gait(t):
    '''return x leg position vector in local coordinates'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    supporting = t_phi <= support_ratio
    x = np.zeros(4)
    x_stance = stride/2 - stride*(t_phi/support_ratio)
    x_swing = -stride/2 + stride*((t_phi - support_ratio)/(1. - support_ratio))
    x = x_swing
    x[supporting] = x_stance[supporting]
    x[t<=T_c] = 0
    return x


def get_y_gait(t):
    '''return y leg position vector in local coordinates'''
    if (type(t) == int) or (type(t) == float):
        y =  np.zeros(4)
    else:
        y = np.zeros((len(t), 4))
    return y


def get_z_gait(t):
    '''return z leg position vector in local coordinates'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    supporting = t_phi <= support_ratio
    z = np.zeros(t_mat.T.shape)
    z_swing = swing_spline((t_phi - support_ratio)/(1. - support_ratio))
    z[~supporting] = z_swing[~supporting]
    z[t<=T_c] = 0
    return z


def get_legstate_gait(t):
    '''return legstate as bool'''
    t_mat = np.tile(t/T_c, (4,1))   # vectorization
    t_phi = (phase+t_mat.T)%1.
    legstate = t_phi <= support_ratio
    legstate[t<=T_c] = True
    return legstate


def get_optimized_gait(t):
    '''returns 3x4 leg coordinates considering optimized leg position'''
    
    # calculate ZMP
    t_ref = np.linspace(0, t[-1], N_ref)
    x_ref = get_x_gait(t_ref) + np.tile(leg_pos[0], (N_ref, 1))
    y_ref = get_y_gait(t_ref) + np.tile(leg_pos[1], (N_ref, 1))
    leg_state = get_legstate_gait(t_ref)
    p_ref = np.zeros((2, N_ref))
    for n in range(N_ref):
        p_ref[0, n] = np.sum(x_ref[n, leg_state[n]])/np.sum(leg_state[n])
        p_ref[1, n] = np.sum(y_ref[n, leg_state[n]])/np.sum(leg_state[n])
        
    # calculate ZMP spline
    spline_com, spline_zmp = zmp.optimized_COM_spline(
                                     p_ref, T_hor=t[-1], 
                                     alpha_x=1, beta_x=0, gamma_x=0,
                                     alpha_y=1, beta_y=1e-6, gamma_y=0)
    
    # calculate desired leg points
    p_zmp = spline_com(t)
    x = get_x_gait(t) - np.tile(p_zmp[0, np.newaxis].T, 4) - 0.05
    y = get_y_gait(t) + 1.5*np.tile(p_zmp[1, np.newaxis].T, 4) # leg numbering wrong?
    z = get_z_gait(t) + z_op
    
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
