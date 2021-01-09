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
collision_ellipsoid = CSV.File("collsion-map-results-2021-01-08_21-46-27/collision_ellipsoid.pd")


Q = hcat(collision_ellipsoid.Q0,
         collision_ellipsoid.Q1,
         collision_ellipsoid.Q2,
         collision_ellipsoid.Q3,
         collision_ellipsoid.Q4,
         collision_ellipsoid.Q5,
         collision_ellipsoid.Q6,
         )
mu = collision_ellipsoid.mu



n = 7 # dimension of the space
max_deg_uv = 2 # degree of moment relaxation
num_pieces = 10 # number of linear pieces
num_iterations=20 # number of iterations of the heuristic
weight_lenght= .1 # trade off between minimizing length and rank
random_seed = 4 # random seed used to initialize the heuristic
a = [1.00, 0.4, 0.0, -1.5, 0.00, 1., -0.0]
b = [-1.00, 0.4, 0.0, -1.5, 0.00, 1., -0.0]
edge_size = 3.5 # edgesize of the bounding box where the trajectory lives

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# Obstacle contraints
# here there is a single obstacle: a disk of center [1-2*t, 0] and radius .5
#
mu = (a .+ b) ./ 2
Q = Q .* 0 + I
moving_disk = [
    (t, x) -> (x .- mu)' * Q * (x .- mu) - 1.
]


# compute optimal piece-wise linear trajectory
opt_trajectory = find_path_using_heuristic(n, moving_disk, edge_size, a, b,
    max_deg_uv, num_pieces, solver,
    weight_lenght,
    num_iterations,
    seed=random_seed)
# opt_trajectory is a function that takes t as input,
# and gives the location of the particle at time t.
Plots.plot(hcat(opt_trajectory.(0:.1:1)...)')

A = hcat(opt_trajectory.(0:.01:1)...)
A
size(A)
writedlm( "computed_trajectory_sos.csv",  A, ',')
