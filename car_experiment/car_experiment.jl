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
using Dates
using DynamicPolynomials

include("plot_helpers.jl")

num_cars = 4
car_rad = .2
obs_rad = .2
obs_pos = [[-1, 1], [1, 1],
           [0, 0]
           ]
obs_vel = [[2, -2], [-2, -2],
           [0, 0]
           ]

g_obstacles = [
    (t, x, y) -> (x-(p[1]+v[1]*t))^2 + (y-(p[2]+v[2]*t))^2 - (car_rad+obs_rad)^2
    for (p, v) ∈ zip(obs_pos, obs_vel)
]

g_plots = [
    (t, x, y) -> (x-(p[1]+v[1]*t))^2 + (y-(p[2]+v[2]*t))^2 - (obs_rad)^2
    for (p, v) ∈ zip(obs_pos, obs_vel)
]

n = 2*num_cars # dimension of the space
max_deg_uv = 2 # degree of moment relaxation
num_pieces = 10 # number of linear pieces
num_iterations=20 # number of iterations of the heuristic
weight_lenght= .08 # trade off between minimizing length and rank
random_seed = 0 # random seed used to initialize the heuristic
a = [0., -1., 1., 0., 0., 1., -1., 0.]
b = [0., 1., -1., 0., 0., -1., 1., 0.]
# a = a .+ rand(Float32, 8) ./ 5.
edge_size = 1.2 # edgesize of the bounding box where the trajectory lives

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# Obstacle contraints
# here there is a single obstacle: a disk of center [1-2*t, 0] and radius .5

obstacles = [
    [(t, x) -> (x[2*(i-1)+1] - x[2*(j-1)+1])^2 +  (x[2*(i-1)+2] - x[2*(j-1)+2])^2 - (2*car_rad)^2
     for i in 1:num_cars for j in (i+1):num_cars]...,
    [ (t, x) -> g(t, x[2i-1], x[2i])
      for i in 1:num_cars for g ∈ g_obstacles]...,
    # [ (t, x) ->  edge_size^2 - x[i]^2
    # for i in 1:2*num_cars]...
]



# compute optimal piece-wise linear trajectory
@show now()
opt_trajectory = find_path_using_heuristic(n, obstacles, edge_size, a, b,
                                           max_deg_uv, num_pieces, solver,
                                           weight_lenght,
                                           num_iterations,
                                           seed=random_seed)

plot_sequence(opt_trajectory, edge_size, g_plots)
plot_animation(opt_trajectory, edge_size, g_plots)

##############
# RRT
##############

include("rrt.jl")
path_found, node_positions, graph_edges = rrt_find_path_tv(n, a, b, check_collision_tv, check_collision_segment_tv; num_iterations=100000, step_size=.3, radius_goal=2.)
println("@ ", Dates.format(now(), "HH:MM:SS"))
@show path_found
opt_trajectory_rrt = t-> rrt_nodes_to_path_tv(node_positions, graph_edges, t)[1:end-1]
plot_sequence(opt_trajectory_rrt, edge_size, g_plots)


################
# NLP
################
include("knitro.jl")
st, st_opt, opt_trajectory_nlp = shortest_path_nlp(
		n, num_cars, car_rad, obs_pos, obs_vel, obs_rad,
		a, b,
		num_pieces,
		edge_size,
		);

println("@ ", Dates.format(now(), "HH:MM:SS"))
@show st
@show st_opt
plot_sequence(opt_trajectory_nlp, edge_size, g_plots)


writedlm( "csv/computed_trajectory_sos.csv",  hcat(opt_trajectory.(0:.05:1)...), ',')
writedlm( "csv/computed_trajectory_rrt.csv",  hcat(opt_trajectory_rrt.(0:.05:1)...), ',')
writedlm( "csv/computed_trajectory_nlp.csv",  hcat(opt_trajectory_nlp.(0:.05:1)...), ',')
if false
    plot_sequence(opt_trajectory, edge_size, g_plots)
    plot_sequence(opt_trajectory_rrt, edge_size, g_plots)
    plot_sequence(opt_trajectory_nlp, edge_size, g_plots)
end
