"""
Helper functions to plot and in this way validate trajectories
"""

import numpy as np
import matplotlib.pyplot as plt
import trajectories as traj
import zero_moment_point as zmp

N = 30     # number of samples for spline calculation AND plot
T_hor = 3   # horizon length
t = np.linspace(0, T_hor, N)
leg_pos = np.array([[0.2, 0.2, -0.2, -0.2],
                    [0.1, -0.1, 0.1, -0.1],
                    [  0,    0,   0,    0]])

# get leg trajectories
x = traj.get_x_gait(t) + np.tile(leg_pos[0], (N, 1))
y = traj.get_y_gait(t) + np.tile(leg_pos[1], (N, 1))
z = traj.get_z_gait(t) + np.tile(leg_pos[2], (N, 1))
leg_state = traj.get_legstate_gait(t)

# desired COP is mean of leg points
p_ref = np.zeros((2, N))
for n in range(N):
    p_ref[0, n] = np.sum(x[n, leg_state[n]])/np.sum(leg_state[n])
    p_ref[1, n] = np.sum(y[n, leg_state[n]])/np.sum(leg_state[n])
    
# calculate ZMP spline
spline_com, spline_zmp = zmp.optimized_COM_spline(
                                     p_ref, T_hor, 
                                     alpha_x=1, beta_x=0, gamma_x=0,
                                     alpha_y=1, beta_y=0, gamma_y=0)
p_com = spline_com(t)
p_zmp = spline_zmp(t)

### Plot ####################################################################
fig, ax = plt.subplots(3, 1, constrained_layout=True)
fig.suptitle("ZMP optimization")
fig.set_figheight(10)

for n in range(4):
    ax[0].plot(x[:, n], y[:, n], label=(f'$leg {str(n)}$'))
    # plot stance feet:
    #ax[1].plot(t[leg_state[:,n] == 1],
    #         (leg_pos[0,n]*leg_state[:,n])[leg_state[:,n] == 1], 
    #         "o", c="k", alpha=0.25)
    #ax[2].plot(t[leg_state[:,n] == 1],
    #         (leg_pos[1,n]*leg_state[:,n])[leg_state[:,n] == 1], 
    #         "o", c="k", alpha=0.5)
    #ax[1].plot(t, x[n])
    #ax[2].plot(t, y[n])

ax[0].plot(p_ref[0], p_ref[1], "--", label="COMref")
ax[0].plot(p_zmp[0], p_zmp[1], label="ZMP")
ax[0].plot(p_com[0], p_com[1], label="COM")    

ax[1].plot(t, p_ref[0], "--", label="$x_{COMref}$")
ax[1].plot(t, p_zmp[0], label="$x_{ZMP}$")
ax[1].plot(t, p_com[0], label="$x_{COM}$")

ax[2].plot(t, p_ref[1], "--", label="$y_{COM,ref}$")
ax[2].plot(t, p_zmp[1], label="$y_{ZMP}$")
ax[2].plot(t, p_com[1], label="$y_{COM}$")
    
ax[0].legend()
ax[0].axis('equal')
ax[0].set_xlabel("x")
ax[0].set_ylabel("y")
ax[1].legend()
ax[1].set_ylim(-0.4,0.4)
ax[1].set_ylabel("x")
ax[1].set_xlabel("t")
ax[2].legend()
ax[2].set_ylim(-0.2,0.2)
ax[2].set_ylabel("y")
ax[2].set_xlabel("t")
    
