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
	using CSV
	using Dates

	using DataFrames
	using DataFramesMeta
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

Debug mode? $(@bind debug_mode CheckBox(default=true))

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
$(@bind step_size Slider(0.00:.01:2, show_value=true, default=1))


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
	@show obstacles
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

# ╔═╡ e590ec80-f87d-11ea-153d-5d9bed4819a5
∞ = Float64(1e9)

# ╔═╡ 8981d24e-f874-11ea-0685-73fd3c241512
md" # SOS"

# ╔═╡ 50972ccc-f876-11ea-05a5-092c2f6da434
function sos_path_successful(path)
	ts = range(0, stop=1, length=100)
	tol = 1e-3
	LinearAlgebra.norm(path(0.) .- [start_x, start_y]) + LinearAlgebra.norm(path(1.) .- [end_x, end_y]) < tol &&
	
	all( f(t, path(t)) >= 0 for t=ts for f=obstacles)
end


# ╔═╡ fd48399e-f87a-11ea-3e0c-dd911a6373f7
function sos_length_path(path)
	ts = range(0, stop=1, length=100)
	sum( LinearAlgebra.norm(path(ti) - path(tii)) for (ti,tii)=zip(ts[1:end-1], ts[2:end]))
end

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
	md"Size path found $(sos_length_path(opt_trajectory))"
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
begin
	function plot_shortest_path(node_positions, graph_edges)
		g = SimpleGraph(size(node_positions, 1))
		for e=graph_edges
			add_edge!(g, e...)
		end
		
		s = [1]
		t = [size(node_positions, 1)]
		size_path = 0
		path = enumerate_paths(dijkstra_shortest_paths(g, s), t)[1]
		for i=1:size(path,1)-1
			a = node_positions[path[i]]
			b = node_positions[path[i+1]]
			size_path += LinearAlgebra.norm(b - a)
			Plots.plot!( [a[1], b[1]], [a[2], b[2]], lw=3, c=:green)
		end
		size_path
	end
	
	function rrt_length_path(node_positions, graph_edges, path_found)
		if ~ path_found
			return ∞
		end
		g = SimpleGraph(size(node_positions, 1))
		for e=graph_edges
			add_edge!(g, e...)
		end
		
		s = [1]
		t = [size(node_positions, 1)]
		size_path = 0
		path = enumerate_paths(dijkstra_shortest_paths(g, s), t)[1]
		for i=1:size(path,1)-1
			a = node_positions[path[i]]
			b = node_positions[path[i+1]]
			size_path += LinearAlgebra.norm(b - a)
		end
		size_path
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
		
	Length path $(rrt_length_path(node_positions, graph_edges, path_found))
	"""
end

# ╔═╡ be90ecb0-f879-11ea-14dd-edff9635e265
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

# ╔═╡ ced25ee2-f874-11ea-3db4-bfd63ce72a87
md"# Benchmark"

# ╔═╡ 5473ac7a-f8ae-11ea-1a13-41a2ebd0b84f
if debug_mode
		range_num_obs = [10, ]
		range_radius_obs = [.3]# range(0.1, stop=.3, length=2)
		range_obs_seed = [1,]
		range_solver_seed = [1,]
else
	range_num_obs = [5, 10, 15, 20, 25]
	range_radius_obs = [.1, .15, .2, .25, .3]
	range_obs_seed = 1:5
	range_solver_seed = 1:5
end

# ╔═╡ cbaba87a-f876-11ea-0533-756c79f0c51c
begin
	
	
	result_sos = []
	result_rrt = []
	for obs_seed=range_obs_seed
		for solver_seed=range_solver_seed
			for number_obs=range_num_obs
				for radius_obs=range_radius_obs
					current_run = (obs_seed, solver_seed, number_obs, radius_obs)
					@show current_run
			# obstacle definitions
			Random.seed!(obs_seed)
			obs_pos = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
			obstacles = [
				(t, x) -> sum( (x .- obs_pos[i, :]).^2 ) - radius_obs^2 for i=1:number_obs
			]
			
			function sos_path_successful(path, obstacles)
				ts = range(0, stop=1, length=100)
				LinearAlgebra.norm(path(0.) .- [start_x, start_y]) + LinearAlgebra.norm(path(1.) .- [end_x, end_y]) < tol &&
				all( f(t, path(t)) >= 0 for t=ts for f=obstacles)
			end
			
			function check_collision(q)
				all( LinearAlgebra.norm(obs_pos[i, :] .- q) > radius_obs for i=1:number_obs)
			end
			
			function check_collision_segment(q1, q2)
				# check if the segment (q1, q2) intersects with the obstacles
				all( check_collision(α .*q1 .+ (1-α) .* q2) for α=range(0, stop=1, length=100) )
			end
			
			
			# Solve with sos
			Random.seed!(solver_seed)
			local opt_trajectory = find_path_using_heuristic(2, obstacles, world_x, 
				[start_x, start_y], [end_x, end_y],
				deg_relaxation, num_pieces, solver,
				weight_lenght,
				num_iterations,
				seed=solver_seed)
			push!(result_sos, [current_run..., 
							sos_path_successful(opt_trajectory, obstacles),
							sos_length_path(opt_trajectory)])
			
			# Solve with RRT
			Random.seed!(solver_seed)
			path_found, node_positions, graph_edges = rrt_find_path(
				[start_x; start_y], 
				[end_x; end_y],
				check_collision, check_collision_segment;
				num_iterations=num_iterations, step_size=step_size,
				radius_goal= .3)
			push!(result_rrt, [current_run..., path_found,
							rrt_length_path(node_positions, graph_edges, path_found)])
				end
			end
		end
	end
end

# ╔═╡ 8be8ca2e-f87b-11ea-383c-a75157b2cf35
function mean(x)
	sum(x) ./ size(x,1)
end

# ╔═╡ f44e403c-f877-11ea-033d-11e91ca00d4e
begin
	col_names = [:obs_seed, :solver_seed, :num_obs, :radius_obs, :success, :length]
	df_sos = DataFrame(hcat(result_sos...)')
	df_rrt = DataFrame(hcat(result_rrt...)')
	for df=(df_sos, df_rrt)
		rename!(df, names(df) .=> col_names)
	end
	rename!(df_rrt, :success => :success_rrt, :length => :length_rrt)
	rename!(df_sos, :success => :success_sos, :length => :length_sos)
	df_rrt_vs_sos = join(df_rrt, df_sos, on= [:obs_seed, :solver_seed, :num_obs, :radius_obs])

end

# ╔═╡ 1fa3b8fe-f878-11ea-099d-eb018d18af63
comb_rrt_vs_sos = combine(groupby(df_rrt_vs_sos, [:num_obs, :radius_obs, :obs_seed ]), :success_rrt => mean, :length_rrt => minimum,  :success_sos => mean, :length_sos => minimum)

# ╔═╡ 98e7f9a4-f87e-11ea-1d55-f9cee67c2702
begin
	success_rate = combine(groupby(df_rrt_vs_sos, [:num_obs, :radius_obs]), :success_rrt => mean, :success_sos => mean)
	
	function make_heatmap_data(data, x, y, v)
	    xs = unique(data[x])
	    ys = unique(data[y])
	    n = length(xs)
	    m = length(ys)
	    A = zeros((n, m))
	    D1 = Dict(x => i for (i,x) in enumerate(xs))
	    D2 = Dict(x => i for (i,x) in enumerate(ys))
	    for i in 1:size(data, 1)
	        xi = data[i, x]
	        yi = data[i, y]
	        vi = data[i, v]
	        A[D1[xi], D2[yi]] = vi
	    end
	    (xs, ys, A)
	end
	
	heat_plots = []
	
	let (x, y, A) = make_heatmap_data(success_rate, :num_obs, :radius_obs, :success_rrt_mean)
	    heat_rrt = heatmap(string.(x), string.(y), A, 
			seriescolor = :blues, clim=(0,1), xlabel="number of obstacles", ylabel="radius" )
		title!("RRT success rate")
		push!(heat_plots, heat_rrt)
	end
	
	let (x, y, A) = make_heatmap_data(success_rate, :num_obs, :radius_obs, :success_sos_mean)
	    heat_sos = heatmap(string.(x), string.(y), A, 
			seriescolor = :blues, clim=(0,1) )
		title!("SOS success rate")
		push!(heat_plots, heat_sos)
	end
	Plots.plot(heat_plots...)
end

# ╔═╡ 0ceeb71c-f892-11ea-17e1-cf1b167c82a7
begin
	groupby_length = combine(groupby(df_rrt_vs_sos, [:obs_seed, :num_obs, :radius_obs]), :length_sos=>minimum, :length_rrt=>minimum)
	#Plots.scatter(groupby_length[:length_sos_minimum], groupby_length[:length_rrt_minimum], ylims=(2,4), xlims=(2, 4), xlabel="sos", ylabel="rrt", legend=false)
	#Plots.plot!(0:10, 0:10, ls=:dash, color=:red)
	Plots.plot(histogram(groupby_length[:length_sos_minimum], bins=2.5:.4:4.5, label="sos"),
	histogram(groupby_length[:length_rrt_minimum], bins=2.5:.4:4.5, label="rrt"))
end

# ╔═╡ 2b081b48-f893-11ea-0f47-f90895af0ac9
groupby_length

# ╔═╡ 081fd8ba-f883-11ea-2dd2-d3d39e7fcf8b
md"# Save to file"

# ╔═╡ 9f54ed10-f881-11ea-0db8-23fa207bd876
begin
	CSV.write("results_sos_vs_rrt-$(Dates.now()).csv", df_rrt_vs_sos)
end

# ╔═╡ 0a67abf2-f892-11ea-294f-fb4916d40a49
df_rrt_vs_sos

# ╔═╡ Cell order:
# ╟─c3d41f9e-f7de-11ea-2a56-b76588d6ef66
# ╟─79575ad0-f7de-11ea-2e97-97d0955e0194
# ╟─8560dd74-f7de-11ea-3c51-7f36d640d43b
# ╟─90161a4a-f7de-11ea-186d-972892ea3c26
# ╟─504017e4-f875-11ea-2928-a1c877aa64b8
# ╠═96523c9a-f7de-11ea-1dd4-67eaad6f968d
# ╟─dc02ed20-f874-11ea-1b15-15f89bb1afc3
# ╟─e590ec80-f87d-11ea-153d-5d9bed4819a5
# ╟─8981d24e-f874-11ea-0685-73fd3c241512
# ╠═50972ccc-f876-11ea-05a5-092c2f6da434
# ╠═fd48399e-f87a-11ea-3e0c-dd911a6373f7
# ╟─bd7b97e4-f7de-11ea-096f-27a885a176c7
# ╟─5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
# ╟─df5b5110-f7de-11ea-3a48-f15db9b1d873
# ╟─725babbe-f874-11ea-2552-9d3968d8b8d9
# ╟─72811ea2-f875-11ea-2c97-b1b5aa72c872
# ╟─771b7552-f875-11ea-3b2c-d7bb19dac4f1
# ╟─68ed372c-f875-11ea-1e3d-07d9438dc9d5
# ╟─937f9bac-f874-11ea-1d59-77ff8d0aaef0
# ╟─be90ecb0-f879-11ea-14dd-edff9635e265
# ╟─ced25ee2-f874-11ea-3db4-bfd63ce72a87
# ╠═5473ac7a-f8ae-11ea-1a13-41a2ebd0b84f
# ╠═cbaba87a-f876-11ea-0533-756c79f0c51c
# ╟─8be8ca2e-f87b-11ea-383c-a75157b2cf35
# ╠═f44e403c-f877-11ea-033d-11e91ca00d4e
# ╟─1fa3b8fe-f878-11ea-099d-eb018d18af63
# ╟─98e7f9a4-f87e-11ea-1d55-f9cee67c2702
# ╠═0ceeb71c-f892-11ea-17e1-cf1b167c82a7
# ╠═2b081b48-f893-11ea-0f47-f90895af0ac9
# ╟─081fd8ba-f883-11ea-2dd2-d3d39e7fcf8b
# ╟─9f54ed10-f881-11ea-0db8-23fa207bd876
# ╠═0a67abf2-f892-11ea-294f-fb4916d40a49
