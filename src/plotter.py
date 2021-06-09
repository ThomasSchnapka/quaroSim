#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions to plot and in this way validate trajectories
"""

import numpy as np
import matplotlib.pyplot as plt
import trajectories as traj
import zero_moment_point as zmp

N = 20
T_hor = 3
t = np.linspace(0, T_hor, N)
leg_pos = np.array([[0.2, 0.2, -0.2, -0.2],
                    [0.1, -0.1, 0.1, -0.1],
                    [0, 0, 0, 0]])

# get leg trajectories
x = traj.get_x_trot(t) + np.tile(leg_pos[0], (N, 1)).T
y = traj.get_y_trot(t) + np.tile(leg_pos[1], (N, 1)).T
z = traj.get_z_trot(t) + np.tile(leg_pos[2], (N, 1)).T
leg_state = traj.get_legstate_trot(t)

# desired COP is mean of leg points
p_ref = np.zeros((2, N))
for n in range(N):
    p_ref[0, n] = np.sum(x[leg_state[n], n])/np.sum(leg_state[n])
    p_ref[1, n] = np.sum(y[leg_state[n], n])/np.sum(leg_state[n])
    
# calculate ZMP spline
spline_zmp = zmp.optimize_ZMP_spline(p_ref, T_hor, 
                                     alpha=1e4, beta=1e-4, gamma=1e-6)
p_zmp = spline_zmp(t)

### Plot ####################################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.set_figheight(10)

for n in range(4):
    ax1.plot(x[n], y[n], label=('x'+str(n)))
    #ax2.plot(t, x[n])
    #ax3.plot(t, y[n])

ax1.plot(p_ref[0], p_ref[1], "--", label="COMref")
ax1.plot(p_zmp[0], p_zmp[1], label="ZMP")   

ax2.plot(t, p_ref[0], "--", label="x_COMref")
ax2.plot(t, p_zmp[0], label="x_ZMP")
#ax2.plot(t, p_zmp[0]-p_ref[0], "--", label="del")


ax3.plot(t, p_ref[1], "--", label="y_COMref")
ax3.plot(t, p_zmp[1], label="y_ZMP")
#ax3.plot(t, p_zmp[1]-p_ref[1], "--", label="del")
    
ax1.legend()
ax1.axis('equal')
ax2.legend()
ax2.set_ylim(-0.4,0.4)
ax3.legend()
ax3.set_ylim(-0.4,0.4)
    
