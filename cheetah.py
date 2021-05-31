import pybullet as p
import pybullet_data as pd
import time
import numpy as np
from  scipy.interpolate import CubicHermiteSpline as hspline

# array that maps quaro convention indices to MiniCheetah joint numbering
jointNumberConversion = np.array(
    [[5, 1, 13, 9],
     [6, 2, 14, 10],
     [4, 0, 12, 8]])

# spline generation
traj = hspline([0,  2, 2.001, 2.2], 
               [0, 20,   0,   0], 
               [0,  0,   0,   0])


def set_leg_angles(angles):
    '''
    changes the robot joint angles according to quaro convention

    Parameters
    ----------
    angles : 3x4 np.ndarray with angles in DEG

    '''
    angles *= 2.0*np.pi/360.0
    p.setJointMotorControlArray(
        robot,
        jointNumberConversion.flatten(),
        p.POSITION_CONTROL, 
        angles.flatten(),
        forces=9999.*np.ones(12))
    

# set up simulation
try:
    p.resetSimulation()
except:
    pass
p.connect(p.GUI)
p.setGravity(0,0,-9.8)
p.setAdditionalSearchPath(pd.getDataPath())
floor = p.loadURDF("plane.urdf")
robot = p.loadURDF("mini_cheetah_custom/mini_cheetah_no_inertia.urdf", [0,0,0.5])
numJoints = p.getNumJoints(robot)
p.changeVisualShape(robot,-1,rgbaColor=[1,1,1,1])
for j in range(numJoints):
    p.changeVisualShape(robot,j,rgbaColor=[1,1,1,1])
    
dt = 1./250.
p.setTimeStep(dt)
tsteps = np.arange(0, 10, dt)
print("loop start")
for t in tsteps:
    # point camera onto robot
    robPos, _ = p.getBasePositionAndOrientation(robot)
    p.resetDebugVisualizerCamera(cameraDistance=1,
                                 cameraYaw=50, cameraPitch=-35, 
                                 cameraTargetPosition=robPos)
    
    # calculate and set position
    # pos = 30*np.sin(2.0*np.pi*t/1)
    pos = traj(t)
    angles = np.ones((3,4))
    angles[0] *= pos
    angles[1] *= -2.*pos
    set_leg_angles(angles)
    
    # simulation step
    p.stepSimulation()
    time.sleep(dt)
    
print("simulation finished!")
        
        
        


