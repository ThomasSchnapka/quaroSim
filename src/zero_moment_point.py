#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Sources:
    
https://www.youtube.com/watch?v=GnFaLl7qwco
http://myweb.sabanciuniv.edu/selimozel/files/2012/10/Zero-Moment-Point-Based-Pace-Reference-Generation-for-Quadruped-Robots-via-Preview-Control.pdf
"""


import numpy as np
import quadprog

# constants
m = 1.
g = 9.81
z_c = 40.    # COM hight
N = 20       # horizon steps
T_tot = 1.5 # horizon length

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

C = np.array([1, 0, -z_c/g])

# state space matrices extended to number of steps
A_e = np.kron(np.eye(N), A)
B_e = np.kron(np.eye(N), B)
C_e = np.kron(np.eye(N), C)

# state space matrices extended to fit 'extended' optvar containing states AND control input
O_x  = np.zeros((n_x, n_x))
O_xu = np.zeros((n_x, n_u))
O_ux = np.zeros((n_u, n_x))
O_u  = np.zeros((n_u, n_u))
O_p = np.zeros((N*n_p, n_u))
O_c = np.zeros(C_e.shape)

A_q = np.block([[A_e,  O_x,  O_xu, O_xu],
                [O_x,  A_e,  O_xu, O_xu],
                [O_ux, O_ux, O_u,  O_u],
                [O_ux, O_ux, O_u,  O_u]])

B_q = np.block([[O_x,  O_x,  B_e,  O_xu],
                [O_x,  O_x,  O_xu, B_e],
                [O_ux, O_ux, O_u,  O_u],
                [O_ux, O_ux, O_u,  O_u]])

C_q = np.block([[C_e, O_c, O_p, O_p],
                [O_c, C_e, O_p, O_p]])


# weighting factors
S_q = np.eye(2*N)
#S_q[2*n_x:] = 0              # penalize only state values
R = np.eye(2*(n_x + n_u))
#R[:2*n_x] = 0                # penalize only control values


def get_zmp(p):
    # for tests:
    p_ref = np.copy(p.flatten())
    
    # convert matrices into a QP
    G_tild = 2*(C_q.T@S_q@C_q + R)
    a_tild = 2*p_ref.T@S_q@C_q
    I = np.ones(2*(n_x + n_u))
    O = np.zeros((2*(n_x + n_u), 1))
    C_tild = np.vstack( (I-A_q-B_q,
                       -(I-A_q-B_q))).T
    b_tild = np.vstack((O, O))[:, 0]        # index necessary for right shape
    
    # solve QP
    sol = quadprog.solve_qp(G_tild, a_tild, C_tild, b_tild)[0]
    result = C_q@sol
    
    p_zmp = np.zeros(p.shape)
    p_zmp[0] = result[:N]
    p_zmp[1] = result[N:]
    return p_zmp



# for tests
if __name__=="__main__":
    p = np.array([[ 0.625     , -5.        ,  6.25      ,  0.625     , -5.        ],
       [ 0.        ,  3.33333333, -3.33333333,  0.        ,  3.33333333]])
    #print(get_zmp(np.random.rand(2, N)))
    print(get_zmp(p/10))


