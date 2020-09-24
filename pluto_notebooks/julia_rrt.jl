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

# ╔═╡ ad7e4182-f312-11ea-08ab-659aec318939
begin
	using Pkg
	Pkg.activate("..")
	using PlutoUI
	using Plots
	using Random
	using LinearAlgebra
	using LightGraphs
	md"""A bunch of imports"""
end

# ╔═╡ 339b48de-f30c-11ea-2de1-ebe2ba2fe62c
plot_disk(x, y, r, c=:black; num=100) = begin
	 # the small circle
	th = range(0, stop=2π, length=num)
	xc = x .+ r .*cos.(th)
	yc = y .+ r .*sin.(th)
	Plots.plot!(xc,yc,c=c)
end

# ╔═╡ d6b26792-f399-11ea-1d61-8ba4a7dd1d7d
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

# ╔═╡ 63429c02-f309-11ea-3830-21f2dff71fcf

md"""
# Data of the problem

world size:

$(@bind world_x Slider(0.00:.01:1, show_value=true, default=1.)) $(@bind world_y Slider(0.00:.01:1, show_value=true, default=1.))

start: 

$(@bind start_x Slider(-1:.01:1, show_value=true, default=-.99)) $(@bind start_y Slider(-1:.01:1, show_value=true, default=-.98))

goal:

$(@bind end_x Slider(0.00:.01:1, show_value=true, default=1.)) $(@bind end_y Slider(0.00:.01:1, show_value=true, default=1.))

number iterations:
$(@bind num_iterations Slider(1:1000, show_value=true, default=100))


Number obstacles:
$(@bind number_obs Slider(1:100, show_value=true, default=10))


Radius of obstacles:
$(@bind radius_obs Slider(0.00:.01:1, show_value=true, default=.1))


obstacle seed:
$(@bind obs_seed Slider(1:100, show_value=true, default=0))


RRT step size
$(@bind step_size Slider(0.00:.01:2, show_value=true, default=.1))


"""

# ╔═╡ 13137234-f30d-11ea-3fb0-c7fb03f32356
begin
	Random.seed!(obs_seed)
	obs_pos = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	
	
	check_collision = q -> all( LinearAlgebra.norm(obs_pos[i, :] .- q) > radius_obs for i=1:number_obs)
	
	# check if the segment (q1, q2) intersects with the obstacles
	check_collision_segment = (q1, q2) -> all( check_collision(α .*q1 .+ (1-α) .* q2) for α=range(0, stop=1, length=100) )
		
	md"obstacles ($number_obs)  generated of radius $(radius_obs)"
end

# ╔═╡ fdea50d0-f30c-11ea-1a0e-a3692e19f287
function rrt_find_path(start, goal, check_collision; num_iterations=10, step_size=.2, radius_goal=.1)
	
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

# ╔═╡ f8672662-f30f-11ea-2380-8ba1806c0a84
begin
	path_found, node_positions, graph_edges = rrt_find_path([start_x; start_y], 
			[end_x; end_y],
			check_collision, num_iterations=num_iterations, step_size=step_size,
		radius_goal= .3)
	md"""RRT finished!!!
	
	Path found? $path_found"""
end

# ╔═╡ 10ad4972-f310-11ea-248d-79a4fdcadc5a
begin
	tol = 1e-2
	plot_path = plot(xlim=(-world_x-tol, world_x+tol), ylim=(-world_x-tol, world_y+tol), 
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

# ╔═╡ 938d283c-f8ef-11ea-2e20-49f7ed6431c5
a = 10

# ╔═╡ Cell order:
# ╟─ad7e4182-f312-11ea-08ab-659aec318939
# ╟─339b48de-f30c-11ea-2de1-ebe2ba2fe62c
# ╟─fdea50d0-f30c-11ea-1a0e-a3692e19f287
# ╟─d6b26792-f399-11ea-1d61-8ba4a7dd1d7d
# ╟─63429c02-f309-11ea-3830-21f2dff71fcf
# ╟─13137234-f30d-11ea-3fb0-c7fb03f32356
# ╟─f8672662-f30f-11ea-2380-8ba1806c0a84
# ╠═10ad4972-f310-11ea-248d-79a4fdcadc5a
# ╟─938d283c-f8ef-11ea-2e20-49f7ed6431c5
