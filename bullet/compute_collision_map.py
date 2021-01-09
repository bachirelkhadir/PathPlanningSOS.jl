"""
Discretize joint space and compute collision map
"""

import numpy as np
import pandas as pd
import os, datetime
from tqdm import tqdm
import pybullet as p
import time
import pybullet_data
from shutil import copyfile

RESULT_DIR = os.path.join(os.getcwd(), "collsion-map-results-"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
RESULT_DIR
os.makedirs(RESULT_DIR)
copyfile("compute_collision_map.py", os.path.join(RESULT_DIR, "compute_collision_map.py"))
csv_file = os.path.join(RESULT_DIR, "collision_map.pd")
print("Result dir: ", RESULT_DIR)

 # #lower limits for null space
 #    self.ll = [-.967, -2, -2.96, 0.19, -2.96, -2.09, -3.05]
 #    #upper limits for null space
 #    self.ul = [.967, 2, 2.96, 2.29, 2.96, 2.09, 3.05]

def cartesian_product(*arrays):
    ndim = len(arrays)
    return np.stack(np.meshgrid(*arrays), axis=-1).reshape(-1, ndim)

NUM_DISC_POINTS = 10
RANGE_JOINTS = [
    (-1, 1, NUM_DISC_POINTS),
    (-2, 2, NUM_DISC_POINTS),
    (-2.96, 2.96, NUM_DISC_POINTS),
    (-0.19, 0.19, NUM_DISC_POINTS),
    (-2.96, 2.96, NUM_DISC_POINTS),
    (-2.09, 2.09, NUM_DISC_POINTS),
    (-3.05, 3.05, NUM_DISC_POINTS),
]

START_POS = [1.00, 0.4, 0.0, -1.5, 0.00, 1., -0.0]
END_POS = [-1.00, 0.4, 0.0, -1.5, 0.00, 1., -0.0]


OBSTACLE_POSE = ([.5, 0, 0.5], p.getQuaternionFromEuler([90, 0.,0]))

OBSTACLE_SCALE = 3

GRID_JOINTS = cartesian_product(*list(map(lambda u: np.linspace(*u), RANGE_JOINTS)))
grid_joints_collision = np.array([False] * len(GRID_JOINTS))


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

p.setRealTimeSimulation(0) # move simulation manually
#############################

##################################
# Collision detection
def dist_robot_obstacle():
    d = 100
    for link in range(7):
        for c in p.getClosestPoints(robot, obstacle_id, distance=100, linkIndexA=link):
            pos_on_A, pos_on_B = c[5:7]
            dist = c[8]
            # p.addUserDebugLine(pos_on_A, pos_on_B, [25, 255, 0], lineWidth=10, lifeTime=1)
            d = min(d, dist)
    return d


def test_collision():
    return  dist_robot_obstacle() < 1e-2


def debug_collision():
    while True:
        time.sleep(1/240.)
        if test_collision():
            p.addUserDebugLine([2, 2, 2], [0, 0, 0], [25, 255, 0], lineWidth=10, lifeTime=1)
        p.stepSimulation()

def test_collision_configuration(joint_angles, pos_obstacle):
    for i, j in enumerate(joint_angles):
        p.resetJointState(robot, i, j)
    p.resetBasePositionAndOrientation(obstacle_id, pos_obstacle, [0, 0, 0, 1])
    return test_collision()


for i, j in tqdm(enumerate(GRID_JOINTS), total=len(GRID_JOINTS)):
    grid_joints_collision[i] = test_collision_configuration(j, OBSTACLE_POSE[0])

collision_map =  pd.DataFrame(GRID_JOINTS, columns=["j{i}".format(i=i) for i in range(7)])
collision_map['collide'] = grid_joints_collision
print("Saving to: ", csv_file) 
collision_map.to_csv(csv_file)
exit()



test_collision_configuration(START_POS, OBSTACLE_POSE[0])
test_collision_configuration([0., 0.4, 0.0, -1.5, 0.0, 1.0, -0.0], OBSTACLE_POSE[0])
test_collision_configuration(END_POS, OBSTACLE_POSE[0])

##################################
p.resetJointState(robot, 1, 0.)
p.resetBasePositionAndOrientation(obstacle_id, [0, 1, 1.5], [0, 0, 0, 1])
dist_robot_obstacle(), test_collision()

p.getClosestPoints(robot, robot, distance=100, linkIndexA=0, linkIndexB=1)[0][8]

collided = p.addUserDebugParameter('Min distance', 0., 1, min(min_dist()))

p.addUserDebugLine([0, 0, 0], [1, 2, 1], [25, 255, 0], lineWidth=10, lifeTime=1)
