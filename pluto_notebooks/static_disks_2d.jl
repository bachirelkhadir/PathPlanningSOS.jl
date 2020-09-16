### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 79575ad0-f7de-11ea-2e97-97d0955e0194
begin
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
	md"""A bunch of imports"""
end

# ╔═╡ c3d41f9e-f7de-11ea-2a56-b76588d6ef66
md"# Disks moving with constant velocity"

# ╔═╡ 8560dd74-f7de-11ea-3c51-7f36d640d43b
md"""
# Data of the problem

world size:

$(@bind world_x Slider(0.00:.01:1, show_value=true, default=1.)) $(@bind world_y Slider(0.00:.01:1, show_value=true, default=1.))


start: 

$(@bind start_x Slider(-1:.01:1, show_value=true, default=-.99)) $(@bind start_y Slider(-1:.01:1, show_value=true, default=-.98))

goal:

$(@bind end_x Slider(0.00:.01:1, show_value=true, default=1.)) $(@bind end_y Slider(0.00:.01:1, show_value=true, default=1.))

number iterations:
$(@bind num_iterations Slider(1:100, show_value=true, default=20))

number pieces:
$(@bind num_pieces Slider(1:100, show_value=true, default=5))


Trade off between minimizing length and rank
$(@bind weight_lenght Slider(0.00:.01:1, show_value=true, default=.1))


Number obstacles:
$(@bind number_obs Slider(1:100, show_value=true, default=10))


Radius of obstacles:
$(@bind radius_obs Slider(0.00:.01:1, show_value=true, default=.1))


obstacle seed:
$(@bind obs_seed Slider(1:100, show_value=true, default=0))


solver seed:
$(@bind solver_seed Slider(1:100, show_value=true, default=0))

"""

# ╔═╡ 90161a4a-f7de-11ea-186d-972892ea3c26
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# ╔═╡ 96523c9a-f7de-11ea-1dd4-67eaad6f968d
begin
	Random.seed!(obs_seed)
	obs_pos = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obs_vel = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obstacles = [
    	(t, x) -> sum( (x .- obs_pos[i, :] .+ t .* obs_vel[i, :]).^2 ) - radius_obs^2 for i=1:number_obs
	]
	
	md"Obstacle definitions here"
end

# ╔═╡ bd7b97e4-f7de-11ea-096f-27a885a176c7
begin
	# compute optimal piece-wise linear trajectory
	Random.seed!(solver_seed)
	opt_trajectory = find_path_using_heuristic(2, obstacles, world_x, 
		[start_x, start_y], [end_x, end_y],
	    2, num_pieces, solver,
	    weight_lenght,
	    num_iterations,
	    seed=solver_seed)
end

# ╔═╡ 5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
md"t $(@bind plot_at_time Slider(0:.01:1, show_value=true, default=0.))"

# ╔═╡ df5b5110-f7de-11ea-3a48-f15db9b1d873
begin
	q = PyPlot.figure()
	PathPlanningSOS.plot_at_time(plot_at_time, world_x, [start_x, start_y], [end_x, end_y],
						obstacles, opt_trajectory)
	PyPlot.title("t = $plot_at_time")
	PyPlot.axes().set_aspect("equal")
	q
end

# ╔═╡ Cell order:
# ╠═c3d41f9e-f7de-11ea-2a56-b76588d6ef66
# ╠═79575ad0-f7de-11ea-2e97-97d0955e0194
# ╟─8560dd74-f7de-11ea-3c51-7f36d640d43b
# ╠═90161a4a-f7de-11ea-186d-972892ea3c26
# ╟─96523c9a-f7de-11ea-1dd4-67eaad6f968d
# ╠═bd7b97e4-f7de-11ea-096f-27a885a176c7
# ╟─5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
# ╟─df5b5110-f7de-11ea-3a48-f15db9b1d873
