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
	using LightGraphs
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

Degree of relaxation
$(@bind deg_relaxation Slider(2:8, show_value=true, default=2))

Number obstacles:
$(@bind number_obs Slider(1:100, show_value=true, default=10))


Radius of obstacles:
$(@bind radius_obs Slider(0.00:.01:1, show_value=true, default=.1))


obstacle seed:
$(@bind obs_seed Slider(1:100, show_value=true, default=0))


solver seed:
$(@bind solver_seed Slider(1:100, show_value=true, default=0))

RRT step size
$(@bind step_size Slider(0.00:.01:1, show_value=true, default=.1))


"""

# ╔═╡ 90161a4a-f7de-11ea-186d-972892ea3c26
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# ╔═╡ 504017e4-f875-11ea-2928-a1c877aa64b8
md"# Obstacle definitions here"

# ╔═╡ 96523c9a-f7de-11ea-1dd4-67eaad6f968d
begin
	Random.seed!(obs_seed)
	obs_pos = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obs_vel =  2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obs_vel .*= 0.
	obstacles = [
    	(t, x) -> sum( (x .- obs_pos[i, :] .+ t .* obs_vel[i, :]).^2 ) - radius_obs^2 for i=1:number_obs
	]
	

	
	
end

# ╔═╡ dc02ed20-f874-11ea-1b15-15f89bb1afc3
begin
	function check_collision(q)
		all( LinearAlgebra.norm(obs_pos[i, :] .- q) > radius_obs for i=1:number_obs)
	end
	function check_collision_segment(q1, q2)
		# check if the segment (q1, q2) intersects with the obstacles
		all( check_collision(α .*q1 .+ (1-α) .* q2) for α=range(0, stop=1, length=100) )
	end
end

# ╔═╡ 8981d24e-f874-11ea-0685-73fd3c241512
md" # SOS"

# ╔═╡ bd7b97e4-f7de-11ea-096f-27a885a176c7
begin
	# compute optimal piece-wise linear trajectory
	Random.seed!(solver_seed)
	opt_trajectory = find_path_using_heuristic(2, obstacles, world_x, 
		[start_x, start_y], [end_x, end_y],
	    deg_relaxation, num_pieces, solver,
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
	PyPlot.plot(
		[i * world_x for i=(-1,1,1,-1,-1)], 
		[j*world_y for j=(-1,-1,1,1,-1)], 
		ls="--", color="red")
	PyPlot.title("SOS @ t = $plot_at_time")
	PyPlot.axes().set_aspect("equal")
	q
end

# ╔═╡ 725babbe-f874-11ea-2552-9d3968d8b8d9
md"# RRT"

# ╔═╡ 72811ea2-f875-11ea-2c97-b1b5aa72c872
plot_disk(x, y, r, c=:black; num=100) = begin
	 # the small circle
	th = range(0, stop=2π, length=num)
	xc = x .+ r .*cos.(th)
	yc = y .+ r .*sin.(th)
	Plots.plot!(xc,yc,c=c)
end

# ╔═╡ 771b7552-f875-11ea-3b2c-d7bb19dac4f1
function rrt_find_path(start, goal, check_collision, check_collision_segment; num_iterations=10, step_size=.2, radius_goal=.1)
	
	node_positions = [ start ]
	graph_edges = []
	Random.seed!(obs_seed)
	
	for counter=1:num_iterations
		
		# pick a random point
		new_x = 2 .* Random.rand(Float32, 2) .- 1
		if ~ check_collision(new_x)
			continue
		end
			
		# find the nearest point in the tree to that point
		dist_to_nodes = [LinearAlgebra.norm(new_x .- n) for n in node_positions]
		closest_node = argmin(dist_to_nodes)
		closest_node_pos = node_positions[closest_node]

		# take step_size in that direction and add that point to the graph
		new_pos = closest_node_pos + step_size .* (new_x-closest_node_pos)
		if ~ check_collision_segment(closest_node_pos, new_pos)
			continue
		end
		push!(node_positions, new_pos )
		push!(graph_edges, (closest_node, size(node_positions, 1)))

		# if we reached the goal
		if LinearAlgebra.norm(new_pos .- goal) < radius_goal
			push!(node_positions, goal )
			push!(graph_edges, (size(node_positions, 1)-1, size(node_positions, 1)))
			return (true, node_positions, graph_edges)
		end
	end
	return (false, node_positions, graph_edges)
end


# ╔═╡ 68ed372c-f875-11ea-1e3d-07d9438dc9d5
function plot_shortest_path(node_positions, graph_edges)
	g = SimpleGraph(size(node_positions, 1))
	for e=graph_edges
		add_edge!(g, e...)
	end
	
	s = [1]
	t = [size(node_positions, 1)]
	
	path = enumerate_paths(dijkstra_shortest_paths(g, s), t)[1]
	for i=1:size(path,1)-1
		a = node_positions[path[i]]
		b = node_positions[path[i+1]]
		Plots.plot!( [a[1], b[1]], [a[2], b[2]], lw=3, c=:green)
	end
end

# ╔═╡ 937f9bac-f874-11ea-1d59-77ff8d0aaef0
begin
	Random.seed!(solver_seed)
path_found, node_positions, graph_edges = rrt_find_path([start_x; start_y], 
		[end_x; end_y],
		check_collision, check_collision_segment;
		num_iterations=num_iterations, step_size=step_size,
		radius_goal= .3)
	md"""RRT finished!!!
	
	Path found? $path_found"""
end

# ╔═╡ ced25ee2-f874-11ea-3db4-bfd63ce72a87
begin
	tol = 1e-2
	plot_path = Plots.plot(xlim=(-world_x-tol, world_x+tol), ylim=(-world_x-tol, world_y+tol), 
		legend=false, aspect_ratio=1,
	title="RRT algorithm")
	for i=1:number_obs
		plot_disk(obs_pos[i, 1], obs_pos[i, 2], radius_obs, :red)
	end
	plot_disk(start_x, start_y, .01, :green)
	plot_disk(end_x, end_y, .01, :green)
	
	for n=node_positions
		plot_disk(n..., .01, :blue)
	end
	for e=graph_edges
		n1 = node_positions[e[1]]
		n2 = node_positions[e[2]]
		Plots.plot!([n1[1], n2[1]], [n1[2], n2[2]], color=:blue)
	end
	plot_shortest_path(node_positions, graph_edges)
	plot_path
end

# ╔═╡ Cell order:
# ╟─c3d41f9e-f7de-11ea-2a56-b76588d6ef66
# ╠═79575ad0-f7de-11ea-2e97-97d0955e0194
# ╟─8560dd74-f7de-11ea-3c51-7f36d640d43b
# ╟─90161a4a-f7de-11ea-186d-972892ea3c26
# ╟─504017e4-f875-11ea-2928-a1c877aa64b8
# ╟─96523c9a-f7de-11ea-1dd4-67eaad6f968d
# ╟─dc02ed20-f874-11ea-1b15-15f89bb1afc3
# ╟─8981d24e-f874-11ea-0685-73fd3c241512
# ╟─bd7b97e4-f7de-11ea-096f-27a885a176c7
# ╟─5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
# ╟─df5b5110-f7de-11ea-3a48-f15db9b1d873
# ╟─725babbe-f874-11ea-2552-9d3968d8b8d9
# ╟─72811ea2-f875-11ea-2c97-b1b5aa72c872
# ╟─771b7552-f875-11ea-3b2c-d7bb19dac4f1
# ╟─68ed372c-f875-11ea-1e3d-07d9438dc9d5
# ╟─937f9bac-f874-11ea-1d59-77ff8d0aaef0
# ╟─ced25ee2-f874-11ea-3db4-bfd63ce72a87
