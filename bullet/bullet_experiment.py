#!/usr/bin/env python3

import numpy as np
from tqdm import tqdm
import pybullet as p
import time
import pybullet_data

# Setup pybullet
physicsClient = p.connect(p.DIRECT)
p.setAdditionalSearchPath(pybullet_data.getDataPath())
pybullet_data.getDataPath()

########################
# Setup scene
p.resetSimulation()
plane_id = p.loadURDF("plane.urdf")
obstacle_id = p.loadURDF("teddy_vhacd.urdf", [0, -1, 0], baseOrientation=p.getQuaternionFromEuler([90, 0.,0]), globalScaling=8.)

robot = p.loadURDF("kuka_experimental/kuka_lbr_iiwa_support/urdf/lbr_iiwa_14_r820.urdf",
                   [0, 0, 0],
                   useFixedBase = 1)

#############################

p.setRealTimeSimulation(0) # move simulation manually
# p.setGravity(0, 0, -.01)
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


# debug_collision()

#################################
#################################
def test_collision_configuration(joint_angles, pos_obstacle):
    for i, j in enumerate(joint_angles):
        p.resetJointState(robot, i, j)
    p.resetBasePositionAndOrientation(obstacle_id, pos_obstacle, [0, 0, 0, 1])
    return test_collision()


for j0 in tqdm(np.linspace(0, 1, int(1e5))):
    test_collision_configuration([0.] * 7, [0, 0, 1.5])

p.resetJointState(robot, 1, 0.)
p.resetBasePositionAndOrientation(obstacle_id, [0, 1, 1.5], [0, 0, 0, 1])
dist_robot_obstacle(), test_collision()

p.getClosestPoints(robot, robot, distance=100, linkIndexA=0, linkIndexB=1)[0][8]

collided = p.addUserDebugParameter('Min distance', 0., 1, min(min_dist()))

p.addUserDebugLine([0, 0, 0], [1, 2, 1], [25, 255, 0], lineWidth=10, lifeTime=1)
