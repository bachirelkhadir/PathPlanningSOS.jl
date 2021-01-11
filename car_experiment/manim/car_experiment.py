import sys
import os
import itertools as it
import numpy as np
from  manim import *
import colorsvg
from helper_functions import *
import pandas as pd

config.background_color = WHITE
config.max_files_cached = 10000

UNIT = 3
WAIT_TIME = 1./60.

def load_img(filename):
    if filename.endswith(".svg"):
        return load_svg(filename)
    return ImageMobject(filename)


def load_svg(filename):
    svg = SVGMobject(filename)
    colorsvg.color_svg_like_file(svg)
    return VGroup(*svg)

def color_car(car, color):
    for i, c in enumerate(car):
        if i not in [1, 2, 5, 6]:
            c.set_color(color)
    return car


# def load_csv_files():
COLORS = (RED, BLUE, GREEN, YELLOW_E)

def circle_cars(cars):
    circles = [Circle(color=col, radius=1/5 * UNIT, fill_color=col, fill_opacity=.3).move_to(c) for c, col in zip(cars, COLORS) ]
    return VGroup(*circles)


def circle_pedestrians(pedestrians):
    circles = [Circle(color=BLACK, radius=1/5 * UNIT, fill_color=BLACK, fill_opacity=.3).move_to(p) for p in pedestrians]
    return VGroup(*circles)


class CarExperiment(Scene):
    def construct(self):
        self.add(NumberPlane().scale(UNIT/2))
        arrows = []
        circles = []
        pedestrians = [load_img("French-Pedestrian-Silhouette.svg").scale(.3).set_color(BLACK)
                        for i in range(3)]
        cars = [color_car(load_img("car-svgrepo-com.svg")[3:].scale(.8), c)
                for c in COLORS]

        cars = stack_group_text(cars)
        cars[0].move_to([0, -UNIT, 0]).rotate(-PI/2)
        cars[1].move_to([0, UNIT, 0]).rotate(PI/2)
        cars[2].move_to([UNIT, 0, 0])
        cars[3].move_to([-UNIT,0, 0]).flip()

        car_circles = circle_cars(cars)
        self.add(car_circles)
        circles.extend(car_circles)

        pedestrians = stack_group_text(pedestrians)
        pedestrians[1].move_to([-UNIT, UNIT, 0]).rotate(PI/6).flip()
        pedestrians[2].move_to([UNIT, UNIT, 0]).rotate(PI/6)

        ped_circles = circle_pedestrians(pedestrians)
        self.add(ped_circles)
        circles.extend(ped_circles)

        for i, p in enumerate(pedestrians[1:]):
            arr = Arrow(p, Dot(), buff=MED_LARGE_BUFF).set_color(BLACK).shift(UP+(-1)**i*LEFT/2).scale(.5)
            self.add(arr)
            arrows.append(arr)


        self.add(VGroup(*cars))
        self.add(VGroup(*pedestrians ))


        arr = DashedVMobject(CurvedArrow(cars[0].get_center()+UP/3+RIGHT/3, cars[1].get_center()+DOWN/3+RIGHT/6, angle=PI/5).set_color(RED))
        arr.shift(RIGHT/3)
        self.add(arr)
        arrows.append(arr)

        arr = arr.copy().set_color(BLUE).rotate(PI).shift(1.5*LEFT)
        self.add(arr)
        arrows.append(arr)

        arr = arr.copy().set_color(YELLOW_E).rotate(PI/2).shift(RIGHT+.7*DOWN)
        self.add(arr)
        arrows.append(arr)

        arr = arr.copy().set_color(GREEN).rotate(PI).shift(1.5*UP)
        self.add(arr)
        arrows.append(arr)

        self.cars = cars
        self.pedestrians = pedestrians
        self.arrows = arrows
        self.circles = circles

        self.current_caption = Tex("")

        # Hide ped 0 and 2
        self.remove(pedestrians[0])
        self.remove(pedestrians[2])
        # self.remove(circles[5])
        self.remove(circles[6])
        self.remove(circles[4])
        self.remove(arrows[1])

    def caption(self, message):
        self.remove(self.current_caption)
        self.current_caption = Tex(message).to_corner(RIGHT).set_color(BLACK)
        self.add(self.current_caption)



class Trajectories(CarExperiment):
    def construct(self):
        CarExperiment.construct(self)

        self.wait(WAIT_TIME)
        self.add_trajectory()
        self.wait(WAIT_TIME)

    def add_trajectory(self):
        sos = pd.read_csv("../csv/computed_trajectory_sos.csv", header=None).values[:, ::1] * UNIT
        rrt = pd.read_csv("../csv/computed_trajectory_rrt.csv", header=None).values[:, ::1] * UNIT
        self.remove(*self.arrows)

        ped_pos = np.array([[-1, 1, 0], [1, 1, 0]]) * UNIT
        ped_vel = np.array([[1, -1, 0], [-1, -1, 0]]) * 2 * UNIT
        full_traj = Dot()
        for name_traj, traj in zip(("sos", "rrt"), (sos, rrt)):
            for with_full_traj in (True, False):
                self.remove(full_traj)
                if with_full_traj:
                    full_traj = get_traj(traj)
                    self.add(full_traj)
                for t in range(len(sos[0])):

                    self.caption(name_traj + f"{t}/{len(sos[0])}")
                    # move cars
                    for i in range(4):
                        self.cars[i].move_to([traj[2*i, t], traj[2*i+1, t], 0])
                        self.circles[i].move_to([traj[2*i, t], traj[2*i+1, t], 0])

                    # move pedesterians
                    for ped, circ, p, v in zip(self.pedestrians[1:], self.circles[5:], ped_pos, ped_vel):
                        ped.move_to(p + t/len(sos[0]) * v)
                        circ.move_to(p + t/len(sos[0]) * v)

                    self.wait(WAIT_TIME)


def get_traj(traj):
    paths = []
    for i, c in enumerate(COLORS):
        path = VMobject()
        path.set_points_as_corners(list(zip(traj[2*i, :], traj[2*i+1, :], [0] * len(traj[0])))).set_color(c)
        # path = DashedVMobject(path)
        paths.append(path)
    return VGroup(*paths)
