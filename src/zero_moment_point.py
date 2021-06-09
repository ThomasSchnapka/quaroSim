#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Sources:
    
https://www.youtube.com/watch?v=GnFaLl7qwco
http://myweb.sabanciuniv.edu/selimozel/files/2012/10/Zero-Moment-Point-Based-Pace-Reference-Generation-for-Quadruped-Robots-via-Preview-Control.pdf


The optimization is done using OSQP which requires the following problem
formulation:
    
minimize        0.5 x' P x + q' x
subject to      l <= A x <= u

"""


import numpy as np
import osqp
import time
from scipy.interpolate import CubicHermiteSpline as hspline
from scipy import sparse

# constants
g = 9.81            # gravitation

def optimize_ZMP_spline(
        p_ref,
        T_hor,
        alpha=1e4, beta=1e-4, gamma=1e-6,
        z_c = 0.25
        ):
    '''
    calculate optimal ZMP and and resulting spline. Algorithm assumes equally
    distributed reference points over time T_hor! 

    Parameters
    ----------
    p_ref : (2xN) np.ndarray, containing reference points in x/y dimension
    T_hor : float, optimization horizon length
    T_hor : float, optimization horizon length
    alpha : float, weighting factor for (p_com - p_ref)
    beta :  float, weighting factor for control values
    gamma : float, weighting factor for state values
    z_c : float, COM hight

    Returns
    -------
    spline : optimal ZMP trajectory in 2D scipy spline

    '''
    
    # create optimization problem
    N = p_ref.shape[1]
    T = T_hor/N
    n_u = N    # number of controls to optimize per dimension
    n_x = 3*N  # number of states to optimize per dimension
    n_p = 1    # number of output variables per dimension
    
    
    # original state space matrices
    A = np.array([[1, T, (T**2)/2],
                  [0, 1,        T],
                  [0, 0,        1]])
    
    B = np.array([[(T**3)/6],
                  [(T**2)/2],
                  [T]])
    
    C = np.array([1., 0, -z_c/g])
    
    
    # state space matrices extended to state/control vector
    A_v = np.kron(np.eye(N, k=-1), A)
    B_v = np.kron(np.eye(N, k=-1), B)
    C_v = np.kron(np.eye(N), C)
    
    # state space matrices extended to fit 'combined' optimization variable
    # containing states AND control input
    O_x  = np.zeros((n_x, n_x))
    O_u  = np.zeros((n_u, n_u))
    O_p = np.zeros((N*n_p, n_u))
    O_xu = np.zeros((n_x, n_u))
    O_ux = np.zeros((n_u, n_x))
    O_Nx = np.zeros((N, n_x))
    O_Nu = np.zeros((N, n_u))
    I_u = np.eye(n_u)
    
    A_tild = np.block([[A_v,  O_x,  B_v,  O_xu],
                       [O_x,  A_v,  O_xu, B_v],
                       [O_ux, O_ux, I_u,  O_u],
                       [O_ux, O_ux, O_u,  I_u]])
    
    C_tild = np.block([[C_v, O_Nx, O_Nu, O_Nu],
                       [O_Nx, C_v, O_Nu, O_Nu],
                       [O_ux, O_ux, I_u, O_u],
                       [O_ux, O_ux, O_u, I_u]])
    
    
    # extented C matrix to use acceleration information for spline generation
    C_ext_x = np.array([0, 1, 0])
    C_ext_u = np.array([-z_c/g])
    C_e = np.kron(np.eye(N), C)
    C_v_ext_x = np.kron(np.eye(N), C_ext_x)
    C_v_ext_u = np.kron(np.eye(N), C_ext_u)
    O_c = np.zeros(C_v.shape)
    O_p = np.zeros(C_v_ext_u.shape)
    C_q_ext = np.block([[C_e,        O_c,       O_p,       O_p      ],
                         [O_c,       C_e,       O_p,       O_p      ],
                         [C_v_ext_x, O_c,       C_v_ext_u, O_p      ],
                         [O_c,       C_v_ext_x, O_p,       C_v_ext_u]])
    
    
    # penalization terms
    R = np.eye(2*n_x+2*n_u)
    Q = np.eye(2*N+2*n_u)
    R[:2*n_x] *= gamma
    Q[:2*N] *= alpha 
    Q[2*N:] *= beta
    
    # combine reference values to single vector
    p = np.copy(p_ref.flatten())
    p_comb = np.hstack((p, np.zeros(2*n_u)))
    
    # convert matrices into QP
    P_QP = 2*(C_tild.T@Q@C_tild + R)
    q_QP = -2*(p_comb.T@Q@C_tild)
    I = np.eye(2*(n_x + n_u))
    O = np.zeros((2*(n_x + n_u), 1))
    A_QP = I-A_tild
    
    # initialize QP
    P_QP = sparse.csc_matrix(P_QP)
    A_QP = sparse.csc_matrix(A_QP)
    prob = osqp.OSQP()
    prob.setup(P=P_QP, q=q_QP, A=A_QP, l=O, u=O, verbose=False)

    # solve QP
    tstart = time.time()
    sol = prob.solve().x
    print(f"solver time: {(time.time()-tstart):.4f}")
    
    # calculate resulting single ZMP points
    result_orig = C_tild@sol
    p_zmp_orig = np.zeros((2, N))
    p_zmp_orig[0] = result_orig[:N]
    p_zmp_orig[1] = result_orig[N:2*N]
    
    # calculate resulting ZMP spline using acceleration information
    result_ext = C_q_ext@sol
    p_zmp_x = result_ext[:N]
    p_zmp_y = result_ext[N:2*N]
    d_p_zmp_x = result_ext[2*N:3*N]
    d_p_zmp_y = result_ext[3*N:]
    tsteps = np.linspace(0, T_hor, N)
    spline = hspline(tsteps,
                     [p_zmp_x, p_zmp_y], 
                     [d_p_zmp_x, d_p_zmp_y],
                     axis=1, extrapolate='periodic')
    return spline


def sample_ref_traj(t):
    '''example COM trajectory'''
    p = np.zeros((2, len(t)))
    p[0, t>T_tot/4] = 0.08*np.sign(np.sin(2*np.pi*t[t>T_tot/4]/(T_tot/2)))
    #p[0, t>T_tot/4] = 0.04*np.sin(2*np.pi*t[t>T_tot/4]/(T_tot/4))
    p[0, t>3*T_tot/4] *= 0
    return p
    

# for tests
if __name__=="__main__":
    '''create example reference and plot optimal ZMP in one dimension'''
    import matplotlib.pyplot as plt
    
    #paramters
    N = 10              # horizon steps
    T_tot = 4           # horizon length
    
    
    # generate sample reference trajectory
    t_n = np.linspace(0, T_tot, N)
    p_ref = sample_ref_traj(t_n)
    
    
    # calculate spline containing optimal ZMP
    spline = optimize_ZMP_spline(p_ref, T_tot)
    
    
    # plot optimal ZMP
    t = np.linspace(0, T_tot, 100)
    vals = spline(t)
    plt.plot(t, vals[0], label="optimal ZMP")
    
    
    # plot sample reference trajectory
    N_plt = 50
    t_plt = np.linspace(0, T_tot, N_plt)
    p_ref = sample_ref_traj(t_plt)
    plt.plot(t_plt, p_ref[0], '--', label="reference")
    plt.legend()


