using Pkg
Pkg.activate("..")
using Revise
using PathPlanningSOS
using MosekTools
using JuMP
using PlutoUI
using Plots
using Random
using PyPlot
using LinearAlgebra
using DataFrames
using CSV
using DelimitedFiles
using ProgressMeter
num_cars = 4
radius_car = .2
rad_obs = .3
g_obstacle = (t, x, y) -> (x-t)^2 + (y+t)^2 - rad_obs^2
n = 2*num_cars # dimension of the space
max_deg_uv = 2 # degree of moment relaxation
num_pieces = 10 # number of linear pieces
num_iterations=20 # number of iterations of the heuristic
weight_lenght= .1 # trade off between minimizing length and rank
random_seed = 0 # random seed used to initialize the heuristic
a = [0., -1., 1., 0., 0., 1., -1., 0.]
b = [0., 1., -1., 0., 0., -1., 1., 0.]
# a = a .+ rand(Float32, 8) ./ 5.
edge_size = 1.2 # edgesize of the bounding box where the trajectory lives

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# Obstacle contraints
# here there is a single obstacle: a disk of center [1-2*t, 0] and radius .5

obstacles = [
    [(t, x) -> (x[2*(i-1)+1] - x[2*(j-1)+1])^2 +  (x[2*(i-1)+2] - x[2*(j-1)+2])^2 - (2*radius_car)^2
    for i in 1:num_cars for j in (i+1):num_cars]...,
    [ (t, x) -> g_obstacle(t, x[2i-1], x[2i])
    for i in 1:num_cars]...,
    # [ (t, x) ->  edge_size^2 - x[i]^2
    # for i in 1:2*num_cars]...
]



# compute optimal piece-wise linear trajectory
opt_trajectory = find_path_using_heuristic(n, obstacles, edge_size, a, b,
                                           max_deg_uv, num_pieces, solver,
                                           weight_lenght,
                                           num_iterations,
                                           seed=random_seed)

# opt_trajectory is a function that takes t as input,
# and gives the location of the particle at time t.
#
fig = figure("Car movement",figsize=(12,2)) # Create a new blank figure
N = 4
car_colors = ["r", "g", "b", "m"]

full_trajectory = hcat(opt_trajectory.(0:.05:1)...)

for t=0:(N-1)
    subplot(101 + N*10 + t, aspect="equal")
    PyPlot.title("t = $(t/(N-1))")
    PyPlot.xlim(-edge_size*1.1, edge_size*1.1)
    PyPlot.ylim(-edge_size*1.1, edge_size*1.1)
    # PyPlot.axis("equal")
    xt = opt_trajectory.(1/(N-1) * t)

    # plt.gcf().gca().add_artist(plt.Circle((0, 0), radius_obs, fill=true))
    PathPlanningSOS.plot_levelset((x, y) -> g_obstacle(1/(N-1) * t, x, y))
    for i=1:num_cars
        plt.gcf().gca().add_artist(plt.Circle(xt[2*i-1:2*i], radius_car, color=car_colors[i], fill=true))
        if t == 0
            PyPlot.plot(full_trajectory[2*i-1, :], full_trajectory[2*i, :], alpha=.5, color=car_colors[i], ls="--")
        end
    end
end


@info "Generating animation"
@showprogress for (i,t)=enumerate(0:.05:1)
    PyPlot.figure()
    PyPlot.axis("equal")
    PyPlot.title("t = $t")

    PyPlot.xlim(-edge_size*1.1, edge_size*1.1)
    PyPlot.ylim(-edge_size*1.1, edge_size*1.1)
    xt = opt_trajectory.(t)


    PathPlanningSOS.plot_levelset((x, y) -> g_obstacle(t, x, y))
    for i=1:num_cars
        plt.gcf().gca().add_artist(plt.Circle(xt[2*i-1:2*i], radius_car, color=car_colors[i], fill=true))
    end

    fn_index = lpad(i, 3, "0")
    PyPlot.savefig("imgs/path_planning_frame_$(fn_index).png")
end

if false

A = hcat(opt_trajectory.(0:.01:1)...)

A
size(A)
writedlm( "computed_trajectory_sos.csv",  A, ',')

end
