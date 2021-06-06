#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Sources:
    
https://www.youtube.com/watch?v=GnFaLl7qwco
http://myweb.sabanciuniv.edu/selimozel/files/2012/10/Zero-Moment-Point-Based-Pace-Reference-Generation-for-Quadruped-Robots-via-Preview-Control.pdf
"""


import numpy as np
import quadprog
import time
from scipy.interpolate import CubicHermiteSpline as hspline
import cvxpy as cp

# constants
g = 9.81
z_c = 0.25#0.615   # COM hight
N = 10       # horizon steps
T_tot = 4  # horizon length

T = T_tot/N
n_u = N    # total number of controls to optimize per dimension
n_x = 3*N  # total number of states to optimize per dimension
n_p = 1    # total number of output variables per dimension


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

# state space matrices extended to fit 'extended' optvar containing states AND control input
O_x  = np.zeros((n_x, n_x))
O_xu = np.zeros((n_x, n_u))
O_ux = np.zeros((n_u, n_x))
O_u  = np.zeros((n_u, n_u))
#O_p = np.zeros((N*n_p, n_u))
#O_c = np.zeros(C_e.shape)
O_p = np.zeros((N*n_p, n_u))
R_p = np.eye(n_u)
O_R = np.zeros((n_u, n_x))
I_u = np.eye(n_u)


A_tild = np.block([[A_v,  O_x,  B_v,  O_xu],
                   [O_x,  A_v,  O_xu, B_v],
                   [O_ux, O_ux, I_u,  O_u],
                   [O_ux, O_ux, O_u,  I_u]])

O_Nx = np.zeros((N, n_x))
O_Nu = np.zeros((N, n_u))
O_ux = np.zeros((n_u, n_x))

C_tild = np.block([[C_v, O_Nx, O_Nu, O_Nu],
                   [O_Nx, C_v, O_Nu, O_Nu],
                   [O_ux, O_ux, I_u, O_u],
                   [O_ux, O_ux, O_u, I_u]])


# including acceleration
C_test_x = np.array([0, 1, 0])
C_test_u = np.array([-z_c/g])
C_e = np.kron(np.eye(N), C)
C_e_test_x = np.kron(np.eye(N), C_test_x)
C_e_test_u = np.kron(np.eye(N), C_test_u)
O_c = np.zeros(C_e.shape)
O_p = np.zeros(C_e_test_u.shape)

C_q_test = np.block([[C_e, O_c, O_p, O_p],
                     [O_c, C_e, O_p, O_p],
                     [C_e_test_x, O_c, C_e_test_u, O_p],
                     [O_c, C_e_test_x, O_p, C_e_test_u]])


# weighting factors

alpha = 1e3#10
beta = 1e-3
gamma = 1e-6#0
R = np.eye(2*n_x+2*n_u)
Q = np.eye(2*N+2*n_u)
R[:2*n_x] *= gamma
Q[:2*N] *= alpha               # penalize only control values
Q[2*N:] *= beta

def get_zmp_spline(p_ref):
    # for tests:
    p = np.copy(p_ref.flatten())
    p_q = np.hstack((p, np.zeros(2*n_u)))
    
    # convert matrices into a QP
    G_QP = 2*(C_tild.T@Q@C_tild + R)
    a_QP = 2*(p_q.T@Q@C_tild).T
    I = np.eye(2*(n_x + n_u))
    O = np.zeros((2*(n_x + n_u), 1))
    C_QP = np.vstack( (I-A_tild,
                     -(I-A_tild))).T
    b_QP = np.vstack((O, O))[:, 0]        # index necessary for right shape
    
    # solve QP
    tstart = time.time()
    sol = quadprog.solve_qp(G_QP, a_QP, C_QP, b_QP)[0]
    print(f"solver time: {(time.time()-tstart):.4f}")
    result_orig = C_tild@sol
    p_zmp_orig = np.zeros((2, N))
    p_zmp_orig[0] = result_orig[:N]
    p_zmp_orig[1] = result_orig[N:2*N]
    '''
    result_test = C_q_test@sol
    p_zmp_x = result_test[:N]
    p_zmp_y = result_test[N:2*N]
    d_p_zmp_x = result_test[2*N:3*N]
    d_p_zmp_y = result_test[3*N:]
    
    tsteps = np.linspace(0, T_tot, N)
    spline = hspline(tsteps,
                     [p_zmp_x, p_zmp_y], 
                     [d_p_zmp_x, d_p_zmp_y],
                     axis=1, extrapolate='periodic')
    return spline, p_zmp_orig, sol
    '''
    return p_zmp_orig, sol

def get_zmp_spline_cvxpy(p_ref):
    # for tests:
    p = np.copy(p_ref.flatten())
    p_q = np.hstack((p, np.zeros(2*n_u)))
    
    # convert matrices into a QP
    G_QP = 2*(C_tild.T@Q@C_tild + R)
    a_QP = 2*(p_q.T@Q@C_tild).T
    I = np.eye(2*(n_x + n_u))
    O = np.zeros((2*(n_x + n_u), 1))
    C_QP = np.vstack( (I-A_tild,
                     -(I-A_tild))).T
    b_QP = np.vstack((O, O))[:, 0]        # index necessary for right shape
    
    # initialize QP
    x = cp.Variable(2*n_x+2*n_u)
    prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(x, G_QP) - a_QP.T @ x),
                     [C_QP.T @ x <= b_QP])
    #prob = cp.Problem(cp.Minimize((C_tild@x-p_q)@Q@(C_tild@x-p_q).T + x.T@R@x),
    #                 [C_QP.T @ x <= b_QP])

    # solve QP
    tstart = time.time()
    prob.solve()
    print(f"solver time: {(time.time()-tstart):.4f}")
    sol = x.value
    result_orig = C_tild@sol
    p_zmp_orig = np.zeros((2, N))
    p_zmp_orig[0] = result_orig[:N]
    p_zmp_orig[1] = result_orig[N:2*N]
    
    result_test = C_q_test@sol
    p_zmp_x = result_test[:N]
    p_zmp_y = result_test[N:2*N]
    d_p_zmp_x = result_test[2*N:3*N]
    d_p_zmp_y = result_test[3*N:]
    
    tsteps = np.linspace(0, T_tot, N)
    spline = hspline(tsteps,
                     [p_zmp_x, p_zmp_y], 
                     [d_p_zmp_x, d_p_zmp_y],
                     axis=1, extrapolate='periodic')
    return p_zmp_orig, sol, spline

def COM_traj(t):
    p = np.zeros((2, len(t)))
    p[0, t>T_tot/4] = 0.08*np.sign(np.sin(2*np.pi*t[t>T_tot/4]/(T_tot/2)))
    #p[0, t>T_tot/4] = 0.04*np.sin(2*np.pi*t[t>T_tot/4]/(T_tot/4))
    p[0, t>3*T_tot/4] *= 0
    return p
    

# for tests
if __name__=="__main__":
    #p = np.array([[ 0.625     , -5.        ,  6.25      ,  0.625     , -5.        ],
    #              [ 0.        ,  3.33333333, -3.33333333,  0.        ,  3.33333333]])
    #p = np.zeros((2, N))
    import matplotlib.pyplot as plt
    t_n = np.linspace(0, T_tot, N)
    p = COM_traj(t_n)
    #p = np.zeros((2, N))
    #p[0] = 0.04*np.sign(np.cos(2*2*np.pi*t_n/T_tot))
    #p[1] = 0.04**np.sign(np.cos(2*2*np.pi*t_n/T_tot))
    #p[0] = 0.04*np.cos(2*np.pi*t_n/T_tot)
    #p[1] = 0.04*np.cos(2*np.pi*t_n/T_tot)
    
    
    #spline, p_orig, sol = get_zmp_spline(p/100)
    #t = np.linspace(0, T_tot, 100)
    #vals = spline(t)
    #plt.plot(t, vals[0])
    #plt.plot(t, vals[1])
    
    # solve with quadprog:
    #p_orig, sol = get_zmp_spline(p)
    #print(p_orig[0])
    t_orig = np.linspace(0, T_tot, N)
    #plt.plot(t_orig, p_orig[0], 'o')
    #plt.plot(t_orig, p_orig[1], 'o')
    
    # solve with cvxpy:
    p_orig, sol, spline = get_zmp_spline_cvxpy(p)
    plt.plot(t_orig, p_orig[0], 'o')
    
    t = np.linspace(0, T_tot, 100)
    vals = spline(t)
    plt.plot(t, vals[0])
    #plt.plot(t, vals[1])
    
    
    N_plt = 50
    t_plt = np.linspace(0, T_tot, N_plt)
    #p = np.zeros((2, N_plt))
    #p[0] = 0.04*np.sign(np.cos(2*2*np.pi*t_plt/T_tot))
    #p[1] = 0.04**np.sign(np.cos(2*2*np.pi*t_plt/T_tot))
    p = COM_traj(t_plt)
    #p[0] = 0.04*np.cos(2*np.pi*t_plt/T_tot)
    #p[1] = 0.04*np.cos(2*np.pi*t_plt/T_tot)
    
    plt.plot(t_plt, p[0], '--')
    #plt.plot(t_plt, p[1]*1e-5, '--')


