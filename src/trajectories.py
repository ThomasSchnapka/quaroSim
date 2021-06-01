# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 21:10:35 2021

@author: z003p2nh
"""
from scipy.interpolate import CubicHermiteSpline as hspline
import numpy as np

# body constants
z_max = 0.389
z_op = 0.25    # usual operating hight

# spline generation
get_y = hspline([0,  2, 2.001, 2.2], 
               [0.2, 0.389,   0.2,   0.2], 
               [0,  0,   0,   0])

get_z = hspline([0,       2,    2.05,    3,    4], 
                [z_op, 0.25,   z_max, z_op, z_op], 
                [0,       0,       0,    0,    0])
