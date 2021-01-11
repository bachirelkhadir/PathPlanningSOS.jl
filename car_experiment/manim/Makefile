##
# Car Example
#
# @file
# @version 0.1



filename = car_experiment.py
scene = CarExperiment
N = 1
q = l

image:
	manim $(filename)  $(scene) -q$q -s -n $N --leave_progress_bars && sleep 1

video:
	manim $(filename)  $(scene) -p -q$q  -n $N  --leave_progress_bars && sleep 1

all:
	manim $(filename)  Trajectories --high_quality  --save_pngs --leave_progress_bars
play:
	xdg-open /home/bachir/Dev/PathPlanningSOS.jl/car_experiment/manim/media/videos/car_experiment/1080p60/Trajectories.mp4

# end
