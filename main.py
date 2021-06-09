'''
dancing cheetah simulation

link to pybullet doc:
https://docs.google.com/document/d/10sXEhzFRSnvFcl3XxNGhnD4N2SedqwdAvK3dsihxVUA/edit#heading=h.mq73m9o2gcpy
'''

import pybullet as p
import pybullet_data as pd
import time
import numpy as np
from src import trajectories as traj

# settings
create_gif= True   # create GIF of simulation
max_force = 10      # maximal joint force

# timing parameters
#dt = 1./250.
#p.setTimeStep(dt)
#tsteps = np.arange(0, 5, dt)
t_end = 6
N_steps = 750
dt = t_end/N_steps
tsteps = np.linspace(0, t_end, N_steps)

# video parameters
fps = 12
N_v = int(t_end*fps)
width = 400
height = 300


# array that maps quaro convention indices to MiniCheetah joint numbering
jointNumberConversion = np.array(
    [[5, 1, 13, 9],
     [6, 2, 14, 10],
     [4, 0, 12, 8]])


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
        forces=max_force*np.ones(12))
    
    
def inverse_kinematics(coordinates):
        '''
        Vectorized inverse kinematics
        
        Parameters
        ----------
        coordinates : (3,4) numpy.ndarray, absolute coordinates for each leg
        
        Returns
        -------
        angles : (3,4) numpy.ndarray
                 angles for each leg in DEG [[femur], [tibia], [coxa]]
        '''
        # local copy of coordinates variable as Python is pass by reference
        coordinates = np.copy(coordinates)
        
        x = coordinates[0]
        y = coordinates[1]
        z = coordinates[2]
        
        # Mini cheetah geometry
        l1 = 0.209
        l2 = 0.18
        g = 0
        h = 0
        
        # inverse kinematics calculation, definitions can be found in doc
        B = np.sqrt(y**2 + z**2)
        A = np.sqrt(B**2 - g**2)# - self.h
        gamma = np.arctan(-y/z) - np.arcsin(g/B)
        C = np.sqrt(A**2 + x**2)
        C1 = ( l1**2 - l2**2 + C**2)/(2*C)
        C2 = (-l1**2 + l2**2 + C**2)/(2*C)
        alpha1 = -np.arctan(-x/A)
        alpha2 = -np.arccos(np.clip(C1/l1, -1, 1))
        teta = -np.arccos(np.clip(C2/l2, -1, 1))
        alpha = alpha1 + alpha2
        beta = -(teta + alpha2)
        
        angles = np.array([alpha, beta, gamma])
        # convert to deg
        angles *= 360.0/(2.0*np.pi)
        
        # replace NaN with zero
        angles = np.nan_to_num(angles)
        
        return angles
    

# caluclate leg positions for whole simulation (closed loop!)
trotPos = traj.get_optimized_trot(tsteps)


# set up simulation
try:
    p.resetSimulation()
except:
    pass
p.connect(p.GUI)
p.setGravity(0,0,-9.8)
p.setAdditionalSearchPath(pd.getDataPath())
floor = p.loadURDF("plane.urdf")
robot = p.loadURDF("src/urdf_files/mini_cheetah_no_inertia.urdf", [0,0,0.4])

# set robot color
numJoints = p.getNumJoints(robot)
p.changeVisualShape(robot,-1,rgbaColor=[1,1,1,1])
for j in range(numJoints):
    p.changeVisualShape(robot,j,rgbaColor=[1,1,1,1]) 
    
# image variables
images = np.zeros((N_v, height, width, 4))
cnt = 0


print("loop start")
for n in range(N_steps):
    t = tsteps[n]
    # point camera onto robot
    robPos, _ = p.getBasePositionAndOrientation(robot)
    p.resetDebugVisualizerCamera(cameraDistance=1,
                                 cameraYaw=50, cameraPitch=-35, 
                                 cameraTargetPosition=robPos)
    # calculate and set position
    pos = trotPos[n]
    angles = inverse_kinematics(pos)
    set_leg_angles(angles)
    
    # save image every now and then
    if create_gif:
        if(n%(int(N_steps/(t_end*fps))) == 0):
            viewMatrix = p.computeViewMatrix(robPos, [0.5, 0.2, 0], [0, 0, 1])
            images[cnt] = p.getCameraImage(width, height, viewMatrix)[2]
            cnt += 1
            cnt = np.min((cnt, N_v-1)) # prevent out of bounds 
        
    # simulation step
    p.stepSimulation()
    time.sleep(dt)
    
print("simulation finished!")

# create GIF
if create_gif:
    print("creating GIF")
    from PIL import Image
    imgs = [Image.fromarray(img.astype(np.uint8)) for img in images]
    imgs[0].save("log.gif", save_all=True, append_images=imgs[1:], fps=fps, loop=0)
        
        
        


