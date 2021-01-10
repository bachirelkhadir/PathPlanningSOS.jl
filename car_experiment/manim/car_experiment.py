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
COLORS = (RED, BLUE, GREEN, YELLOW)

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
        pedestrians = [load_img("French-Pedestrian-Silhouette.svg").scale(.3).set_color(BLACK)
                        for i in range(3)]
        cars = [color_car(load_img("car-svgrepo-com.svg")[3:].scale(.8), c)
                for c in COLORS]

        cars = stack_group_text(cars)
        cars[0].move_to([0, -UNIT, 0]).rotate(-PI/2)
        cars[1].move_to([0, UNIT, 0]).rotate(PI/2)
        cars[2].move_to([UNIT, 0, 0])
        cars[3].move_to([-UNIT,0, 0]).flip()
        self.add(circle_cars(cars))

        pedestrians = stack_group_text(pedestrians)
        pedestrians[1].move_to([-UNIT, UNIT, 0]).rotate(PI/6).flip()
        pedestrians[2].move_to([UNIT, UNIT, 0]).rotate(PI/6)
        self.add(circle_pedestrians(pedestrians))

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


        # self.add_trajectory()
        self.wait()


    def add_trajectory(self):

        sos =pd.read_csv("../csv/computed_trajectory_sos.csv", header=None).values * UNIT
        rrt=pd.read_csv("../csv/computed_trajectory_rrt.csv", header=None).values * UNIT
        for traj in (sos, rrt):
            path = VMobject()
            path.set_points_as_corners(list(zip(traj[0, :], traj[1, :], [0] * len(sos[0]))))
            path.set_color(BLACK)
            
            self.add(path)
