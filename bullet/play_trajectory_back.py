import numpy as np
from tqdm import tqdm
import pybullet as p
import time
import pandas as pd
import pybullet_data
import matplotlib.pyplot as plt
plt.style.use('dark_background')

START_POS = [1.00, 0.4, 0.0, -1.5, 0.00, 1., -0.0]
END_POS = [-1.00, 0.4, 0.0, -1.5, 0.00, 1., -0.0]

START_POS = np.array(START_POS)
END_POS = np.array(END_POS)

OBSTACLE_POSE = ([.5, 0, 0.5], p.getQuaternionFromEuler([90, 0.,0]))
OBSTACLE_SCALE = 3


OBSTACLE_URDF = "teddy_vhacd.urdf"
KUKA_URDF = "kuka_experimental/kuka_lbr_iiwa_support/urdf/lbr_iiwa_14_r820.urdf"


# Setup pybullet
physicsClient = p.connect(p.GUI)
# physicsClient = p.connect(p.DIRECT)
p.setAdditionalSearchPath(pybullet_data.getDataPath())
pybullet_data.getDataPath()

########################
# Setup scene
p.resetSimulation()
# Load data
plane_id = p.loadURDF("plane.urdf")
obstacle_id = p.loadURDF(OBSTACLE_URDF, OBSTACLE_POSE[0], baseOrientation=OBSTACLE_POSE[1], globalScaling=OBSTACLE_SCALE)
robot = p.loadURDF(KUKA_URDF, [0, 0, 0], useFixedBase=1)

#############################

def set_configuration(joint_angles, pos_obstacle):
    for i, j in enumerate(joint_angles):
        p.resetJointState(robot, i, j)
    p.resetBasePositionAndOrientation(obstacle_id, pos_obstacle, [0, 0, 0, 1])



# trajectory = pd.read_csv("computed_trajectory_sos.csv", index).values
trajectory = list(map(lambda s: list(map(float, s.split(","))), open("computed_trajectory_sos.csv").readlines()))
trajectory = np.array(trajectory)

# plt.plot(trajectory.T)
for j in trajectory.T:
    set_configuration(j, OBSTACLE_POSE[0])
    # print(j)
    time.sleep(.03)
