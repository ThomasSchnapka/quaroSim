"""

Sources:
    
https://www.youtube.com/watch?v=GnFaLl7qwco
http://myweb.sabanciuniv.edu/selimozel/files/2012/10/Zero-Moment-Point-Based-Pace-Reference-Generation-for-Quadruped-Robots-via-Preview-Control.pdf


The optimization is done using OSQP which requires the following problem-
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

def optimized_COM_spline(
        p_ref,
        T_hor,
        alpha_x=1e4, beta_x=1e-4, gamma_x=1e-6,
        alpha_y=1e4, beta_y=1e-4, gamma_y=1e-6,
        z_c = 0.35
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
             spline function returns 2xM spline points given t (1xM) time steps

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
                       [O_ux, O_ux, O_u, O_u],
                       [O_ux, O_ux, O_u, O_u]])
    
    
    # extented C matrix to use acceleration information for spline generation
    C_zmp_x = np.array([0, 1, 0])
    C_zmp_u = np.array([-z_c/g])
    C_e = np.kron(np.eye(N), C)
    C_v_zmp_x = np.kron(np.eye(N), C_zmp_x)
    C_v_zmp_u = np.kron(np.eye(N), C_zmp_u)
    O_c = np.zeros(C_v.shape)
    O_p = np.zeros(C_v_zmp_u.shape)
    C_q_zmp = np.block([[C_e,        O_c,       O_p,       O_p      ],
                         [O_c,       C_e,       O_p,       O_p      ],
                         [C_v_zmp_x, O_c,       C_v_zmp_u, O_p      ],
                         [O_c,       C_v_zmp_x, O_p,       C_v_zmp_u]])
    
    
    # penalization terms
    R = np.eye(2*n_x+2*n_u)
    Q = np.eye(2*N+2*n_u)
    Q[   :1*N] *= alpha_x
    Q[  N:2*N] *= alpha_y
    Q[-2*n_u:n_u] *= beta_x
    Q[  -n_u:   ] *= beta_y
    R[     :1*n_x] *= gamma_x
    R[1*n_x:2*n_x] *= gamma_y
    R[2*n_x:     ] *= 0
    
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
    print("starting solver for ZMP optimization")
    tstart = time.time()
    sol = prob.solve().x
    print(f"solver time: {((time.time()-tstart)):.9f}s")
    
    # calculate resulting single ZMP points
    # result_orig = C_tild@sol
    # p_zmp_orig = np.zeros((2, N))
    # p_zmp_orig[0] = result_orig[:N]
    # p_zmp_orig[1] = result_orig[N:2*N]
    
    # calculate resulting COM spline
    p_com_x   = sol[     :n_x  :3]
    d_p_com_x = sol[    1:n_x  :3]
    p_com_y   = sol[n_x  :2*n_x:3]
    d_p_com_y = sol[n_x+1:2*n_x:3]
    tsteps = np.linspace(0, T_hor, N)
    spline_com = hspline(tsteps,
                         [p_com_x, p_com_y], 
                         [d_p_com_x, d_p_com_y],
                         axis=1, extrapolate='periodic')
    
    # calculate resulting ZMP spline using acceleration information
    result_zmp = C_q_zmp@sol
    p_zmp_x = result_zmp[:N]
    p_zmp_y = result_zmp[N:2*N]
    d_p_zmp_x = result_zmp[2*N:3*N]
    d_p_zmp_y = result_zmp[3*N:]
    tsteps = np.linspace(0, T_hor, N)
    spline_zmp = hspline(tsteps,
                         [p_zmp_x, p_zmp_y], 
                         [d_p_zmp_x, d_p_zmp_y],
                         axis=1, extrapolate='periodic')
    
    return spline_com, spline_zmp




