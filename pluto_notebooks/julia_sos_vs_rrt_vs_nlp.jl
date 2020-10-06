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
	using Statistics

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
	using KNITRO 
	using ProgressMeter
	using DynamicPolynomials
	md"""A bunch of imports"""
end

# ╔═╡ c3d41f9e-f7de-11ea-2a56-b76588d6ef66
md"# Disks moving with constant velocity"

# ╔═╡ 98a3f246-05e8-11eb-201d-a133f7585b7c
md"Debug mode? $(@bind debug_mode CheckBox(default=true))"

# ╔═╡ a39d4a4a-05e7-11eb-2f92-e118170808f8
n = 2

# ╔═╡ abd44b0a-05e7-11eb-07ed-c7c73096570e
edge_size = 1.01

# ╔═╡ b1917388-05e7-11eb-32fb-45c9ace941de
a = [-.99 for i=1:n]

# ╔═╡ ba98a5aa-05e7-11eb-06c9-653431e3ec70
b = [.99 for i=1:n]

# ╔═╡ c793b196-05e7-11eb-3fd5-dd7fd6b48e49
num_sos_iterations = 20

# ╔═╡ e956f92a-05ea-11eb-3849-c9c3cdf5fb47
num_rrt_iterations = 10000

# ╔═╡ cb711d12-05e7-11eb-3056-5fa4b9853578
num_pieces = 5

# ╔═╡ cf3ed704-05e7-11eb-042c-bb774ba8306e
weight_lenght = .1

# ╔═╡ d47719ac-05e7-11eb-35b4-ed887c5c643f
deg_relaxation = 2

# ╔═╡ d86e55c8-05e7-11eb-037d-894ea2b56788
number_obs = 10

# ╔═╡ de847bf6-05e7-11eb-39ca-874f2cdd74fc
radius_obs = .2

# ╔═╡ 587fb63e-066d-11eb-3499-0ff6ff490bb2
md"Static obstacles? $(@bind static_mode CheckBox(default=true))"

# ╔═╡ 1d3b7e12-05e8-11eb-3f76-3b7dc80b948f
speed_obs = static_mode ? 0.  : 1.

# ╔═╡ e856be64-05e7-11eb-2260-77c2b187e2ed
obs_seed = 4

# ╔═╡ ec1ce6ea-05e7-11eb-3816-11b297344745
solver_seed = 1

# ╔═╡ f06b8128-05e7-11eb-23db-5fbe3f7cd461
step_size = 2.

# ╔═╡ 90161a4a-f7de-11ea-186d-972892ea3c26
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# ╔═╡ e590ec80-f87d-11ea-153d-5d9bed4819a5
∞ = Float64(1e9)

# ╔═╡ 1631c244-05f1-11eb-330a-4b2bdc5073a0
obstacle_setup = Dict(:n => n, 
	:number_obs => number_obs, 
	:speed_obs => speed_obs,
	:radius_obs => radius_obs,
	:obs_seed => obs_seed,
	:a => a,
	:b => b)

# ╔═╡ 2011598c-05ec-11eb-1890-850db5315bd5

function generate_obstacle_setup(obstacle_setup)
	# reading necessary config variables
	obs_seed = obstacle_setup[:obs_seed]
	n = obstacle_setup[:n]
	number_obs = obstacle_setup[:number_obs]
	speed_obs = obstacle_setup[:speed_obs]
	radius_obs = obstacle_setup[:radius_obs]

	# generating obstacles
	Random.seed!(obs_seed)
	obs_pos = 2 .* Random.rand(Float64, (number_obs, n)) .- 1
	obs_vel = 2  .* Random.rand(Float64, (number_obs, n)) .- 1
	obs_vel .*= speed_obs
	
	obstacles = [
	(t, x) -> sum( (x .- obs_pos[i, :] .- t .* obs_vel[i, :]).^2 ) - radius_obs^2 for i=1:number_obs
	]

	# helper functions to check collisions between segments and obstacles
# 	function sos_path_successful(path, obstacles)
# 		tol = 1e-2
# 		ts = range(0, stop=1, length=100)
# 		LinearAlgebra.norm(path(0.) .- a) + LinearAlgebra.norm(path(1.) .- b) < tol &&
# 		all( f(t, path(t)) >= 0 for t=ts for f=obstacles)
# 	end

	check_collision_tv(q, t)= all( obs(t, q) >= 0 for obs ∈ obstacles) &&
							  all( abs.(q) .<= edge_size)
	

	# check if the segment (q1, q2) intersects with the obstacles
	check_collision_segment_tv = (q1, q2) -> begin
		t1 = q1[end]
		t2 = q2[end]
		ts = range(0, stop=1, length=100)
		t1 < t2 &&
		all( check_collision_tv( (1-α) .*q1[1:end-1] .+ α .* q2[1:end-1], 
		(1-α) *t1 .+ α * t2) for α ∈ ts)
	end
	check_collision(q) = check_collision_tv(q, 0.)

	check_collision_segment(q1, q2) = 
		check_collision_segment_tv([q1..., 0], [q1..., 1])
	
	Dict(:obstacles => obstacles,  
		:obs_pos => obs_pos,
		:obs_vel => obs_vel,
		:radius_obs => radius_obs,
		:check_collision_tv => check_collision_tv, 
		:check_collision_segment_tv => check_collision_segment_tv, 
		:check_collision => check_collision,
		:check_collision_segment => check_collision_segment,
		)

end

# ╔═╡ 4d5b9cae-0605-11eb-2cf3-576448d0439e
md"t $(@bind plot_at_time Slider(0:.01:1, show_value=true, default=0.))"

# ╔═╡ 691cc1b2-0604-11eb-01a1-e9db3148a92a
function plot_obstacles(obstacle_setup, opt_trajectory=(t->a))
	obstacles = generate_obstacle_setup(obstacle_setup)[:obstacles]
	q = PyPlot.figure()
	PathPlanningSOS.plot_at_time(plot_at_time, edge_size, a, b,
						obstacles, opt_trajectory)
	PyPlot.title("t = $plot_at_time")
	PyPlot.axes().set_aspect("equal")
	q
end

# ╔═╡ c6c0c038-0600-11eb-1c80-b51afcbbab50
md"# Path helpers"

# ╔═╡ cb6f6e86-0600-11eb-33ec-416d9a847dc3
function length_path(path)
        ts = range(0, stop=1, length=100)
        sum( LinearAlgebra.norm(path(ti) - path(tii)) for (ti,tii)=zip(ts[1:end-1], ts[2:end]))
end

# ╔═╡ 8981d24e-f874-11ea-0685-73fd3c241512
md" # SOS"

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
function rrt_find_path(dim, start, goal, check_collision, check_collision_segment; num_iterations=10, step_size=.2, radius_goal=.1)
	
	node_positions = [ start ]
	graph_edges = []
	Random.seed!(obs_seed)
	
	for counter=1:num_iterations
		
		# pick a random point
		new_x = 2 .* Random.rand(Float32, dim) .- 1
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


# ╔═╡ 09836820-05ee-11eb-08df-993638bce61f
function rrt_nodes_to_path(node_positions, graph_edges, t)
	g = LightGraphs.SimpleGraph(size(node_positions, 1))
	for e=graph_edges
		add_edge!(g, e...)
	end
	
	xi = [1]
	xf = [size(node_positions, 1)]
	all_paths = enumerate_paths(dijkstra_shortest_paths(g, xi), xf)

	path = all_paths[1]
	if size(path, 1) == 0
 		return node_positions[1]
 	end
 	s = size(path,1)
	
 	curr_segment = Int(floor(t * s, digits=0))+1
	t1 = (curr_segment-1) / s
	t2 = t1 + 1. /s
	a = node_positions[path[max(1, min(s, curr_segment))]]
	b = node_positions[path[max(1, min(s, curr_segment+1))]]
	
	α = 1 - (t - t1) / (t2 - t1)
	c = α .* a .+ (1-α) .* b
	
end


# ╔═╡ 28bd2dd0-05ed-11eb-3389-4dc1b4c7bfe1
md"# RRT TV"

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
	
	
	function rrt_length_path_tv(node_positions, graph_edges, path_found)
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
			size_path += LinearAlgebra.norm(b[1:end-1] .- a[1:end-1])
		end
		size_path
	end
	
end

# ╔═╡ b696d130-05e9-11eb-374c-1b166d28069d
function rrt_find_path_tv(dim, start, goal, check_collision, check_collision_segment; num_iterations=10, step_size=.2, radius_goal=.1)
	
	node_positions = [ [start..., 0] ]
	graph_edges = []
	Random.seed!(obs_seed)
	
	for counter=1:num_iterations
		
		# pick a random point
		new_x = 2 .* Random.rand(Float32, dim) .- 1
		new_t = Random.rand(Float32, 1)[1]
		if ~ check_collision(new_x, new_t)
			continue
		end
			
		# find the nearest point in the tree to that point
		dist_to_nodes = [LinearAlgebra.norm(new_x .- n[1:end-1]) for n in node_positions]
		closest_node = argmin(dist_to_nodes)
		closest_node_pos = node_positions[closest_node]

		# take step_size in that direction and add that point to the graph
		new_pos = zeros(Float32, dim+1)
		for j ∈ 1:size(start, 1)
			new_pos[j] = closest_node_pos[j] + step_size * (new_x[j]-closest_node_pos[j])
		end

		new_pos[end] = new_t
		
		if ~ check_collision_segment(closest_node_pos, new_pos)
			continue
		end
		push!(node_positions, new_pos )
		push!(graph_edges, (closest_node, size(node_positions, 1)))

		# if we reached the goal
		if LinearAlgebra.norm(new_pos[1:end-1] .- goal) < radius_goal
			push!(node_positions, [goal..., new_pos[end]] )
			push!(graph_edges, (size(node_positions, 1)-1, size(node_positions, 1)))
			return (true, node_positions, graph_edges)
		end
	end
	return (false, node_positions, graph_edges)
end

# ╔═╡ 2bf3de7c-05ed-11eb-2189-9b7d947b3c19
function rrt_nodes_to_path_tv(node_positions, graph_edges, t)
	g = LightGraphs.SimpleGraph(size(node_positions, 1))
	for e=graph_edges
		add_edge!(g, e...)
	end
	
	xi = [1]
	xf = [size(node_positions, 1)]
	
	path = enumerate_paths(dijkstra_shortest_paths(g, xi), xf)[1]
	current_pos_found = false
	for i=1:size(path,1)-1
		a = node_positions[path[i]]
		b = node_positions[path[i+1]]
		t1 = a[end]
		t2 = b[end]
				
		if t1 <= t && t <= t2
			current_pos_found = true
			α = 1 - (t - t1) / (t2 - t1)
			c = α .* a .+ (1-α) .* b
		end
	end
	if ~current_pos_found
		c = node_positions[path[end]]
	end
	c
end


# ╔═╡ 43e4f67a-05e4-11eb-29d3-83cd2faaee75
md"# NLP"

# ╔═╡ 6bce6542-05e5-11eb-1322-81821b08c030
function shortest_path_nlp(
		n, obs_pos, obs_vel, radius_obs, 
		a, b,
		num_pieces,
		edge_size,
		)
	@info "Model"
	
	model = Model(optimizer_with_attributes(KNITRO.Optimizer, "honorbnds" => 1,
											"outlev" => 0, "algorithm" => 4,
											))
	@variable(model, u[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	@variable(model, v[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	
	@info "Constraints"
		lengths = @NLexpression model [i=1:size(v,1)] sqrt(sum(v[i, j]^2 for j=1:n))
		
		for j=size(obs_pos, 1)
			for i=1:num_pieces
				for t=0:.1:1
					obs_pos_t = obs_pos[j, :] .+ t .*  obs_vel[j, :]
					@NLconstraint model  sum( (u[i,k] + t*v[i,k]-obs_pos_t[k] )^2 for k=1:n) >= radius_obs^2
				end
			end
		end
		
		for i=1:num_pieces-1
			@constraint model  u[i,:] .+ v[i,:] .== u[i+1,:]
		end
		
		@constraint model u[1,:] .== a
		@constraint model u[end,:] .+ v[end,:] .== b;
		
		
		for j=(1, -1)
			for i=1:num_pieces
				@constraint model j .* u[i,:] .<= edge_size
			end
		end
		@info "Set Objective"
		@NLobjective model Min sum(l for l=lengths)
		@info "Optimize"
		optimize!(model)
		termination_status(model), objective_value(model)
		st = termination_status(model)
		st_opt = objective_value(model)
	
		@polyvar t
		opt_trajectory_pieces = [value.(u[i, :]) + t .* value.(v[i, :]) for i=1:num_pieces]
		(st, st_opt, 
		PathPlanningSOS.pieces_to_trajectory(opt_trajectory_pieces))

end

# ╔═╡ ced25ee2-f874-11ea-3db4-bfd63ce72a87
md"# Benchmark"

# ╔═╡ 5473ac7a-f8ae-11ea-1a13-41a2ebd0b84f
if debug_mode
		range_n = [2, 3, ]
		range_obs_seed = [1, 2]
		range_speed_obs = [1.]
		range_solver_seed = [1,]
else
	range_n = 2:4
	range_obs_seed = 1:10
	range_speed_obs = [0., 1.]
	range_solver_seed = [1,]
end

# ╔═╡ 8be8ca2e-f87b-11ea-383c-a75157b2cf35
function mean(x)
	sum(x) ./ size(x,1)
end

# ╔═╡ d358181e-0600-11eb-1ad9-4316519b7784
function is_path_valid(path, obstacle_setup; tol=1e-2, 
					num_discretization_points=1000,
					max_jump=.1)
	obstacles = generate_obstacle_setup(obstacle_setup)[:obstacles]
	a = obstacle_setup[:a]
	b = obstacle_setup[:b]
	ts = range(0, stop=1, length=num_discretization_points)
	norm =  LinearAlgebra.norm
	start_and_dest_const = norm(path(0.) .- a) + norm(path(1.) .- b) < tol
	collision_avoidance_const = all( f(t, path(t)) >= 0 for t=ts for f=obstacles)
	continuity_const = all( mean(abs.(path(ti) .- path(tii))) <= max_jump for (ti, tii)=zip(ts[1:end-1], ts[2:end]) )
	
	start_and_dest_const && collision_avoidance_const# && continuity_const
end


# ╔═╡ ebc7c88c-05e6-11eb-13cd-455b20227590
begin
	function run_nlp(obstacle_setup)
	generated_obs = generate_obstacle_setup(obstacle_setup)

		
	#obs_vel = obs_pos .* 0
	opt_trajectory_nlp = shortest_path_nlp(obstacle_setup[:n],
						generated_obs[:obs_pos], 
						generated_obs[:obs_vel], 
						generated_obs[:radius_obs], 
						obstacle_setup[:a], obstacle_setup[:b],
						num_pieces,
						edge_size)[end]
	end
	opt_trajectory_nlp = run_nlp(obstacle_setup)
	is_path_valid(opt_trajectory_nlp, obstacle_setup)
end

# ╔═╡ fce9d6a4-0603-11eb-10c7-27fc653d8c99
if n == 2
	q_rrt_nlp = plot_obstacles(obstacle_setup, opt_trajectory_nlp)
	PyPlot.plot(eachcol(hcat(opt_trajectory_nlp.(0:.01:1)...)')...)
	q_rrt_nlp
end

# ╔═╡ 69813a14-0615-11eb-2398-31b3f1b9c905
function path_smoothness(path; dt=0.01)
	P = hcat(path.(0:dt:1)...)
	sum(std.(eachrow(P[:, 1:end-1] .- P[:, 2:end])))
end

# ╔═╡ bd7b97e4-f7de-11ea-096f-27a885a176c7
begin
	function run_sos(obstacle_setup)
		generated_obs = generate_obstacle_setup(obstacle_setup)
		# compute optimal piece-wise linear trajectory
		Random.seed!(solver_seed)
		opt_trajectory_sos = find_path_using_heuristic(obstacle_setup[:n], 
			generated_obs[:obstacles], 
			edge_size, 
			obstacle_setup[:a], obstacle_setup[:b],
			deg_relaxation, num_pieces, solver,
			weight_lenght,
			num_sos_iterations,
			seed=solver_seed)
		opt_trajectory_sos
	end
	opt_trajectory_sos = run_sos(obstacle_setup)
	is_path_valid(opt_trajectory_sos, obstacle_setup), "smoothness", path_smoothness(opt_trajectory_sos)
end

# ╔═╡ b3b554bc-05e8-11eb-30ee-13748ef82f3b
if n == 2
	q_sos = plot_obstacles(obstacle_setup, opt_trajectory_sos)
	PyPlot.plot(eachcol(hcat(opt_trajectory_sos.(0:.01:1)...)')...)
	q_sos
end

# ╔═╡ 3eda12e4-0616-11eb-3791-03110d0f8432
LinearAlgebra.norm(opt_trajectory_sos(0.) .- opt_trajectory_sos(0.01))

# ╔═╡ 937f9bac-f874-11ea-1d59-77ff8d0aaef0
begin
	function run_rrt(obstacle_setup)
	generated_obs = generate_obstacle_setup(obstacle_setup)

		
	Random.seed!(solver_seed)
	path_found, node_positions, graph_edges = rrt_find_path(obstacle_setup[:n],
			obstacle_setup[:a], obstacle_setup[:b],
		generated_obs[:check_collision], generated_obs[:check_collision_segment];
		num_iterations=num_rrt_iterations, step_size=step_size,
		radius_goal= .3)
	opt_trajectory_rrt = t-> rrt_nodes_to_path(node_positions, graph_edges, t)
	end
	opt_trajectory_rrt = run_rrt(obstacle_setup)
	is_path_valid(opt_trajectory_rrt, obstacle_setup), "smoothness", path_smoothness(opt_trajectory_rrt)
	
	
end

# ╔═╡ 8bf7a05a-0602-11eb-0940-1b290811b69d
if n == 2
	q_rrt = plot_obstacles(obstacle_setup, opt_trajectory_rrt)
	PyPlot.plot(eachcol(hcat(opt_trajectory_rrt.(0:.01:1)...)')...)
	q_rrt
end

# ╔═╡ 3408f502-05ed-11eb-1316-cdd8235769bf
begin
	function run_rrt_tv(obstacle_setup)
	generated_obs = generate_obstacle_setup(obstacle_setup)

		
	Random.seed!(solver_seed)
	path_found, node_positions, graph_edges = rrt_find_path_tv(obstacle_setup[:n],
			obstacle_setup[:a], obstacle_setup[:b],
		generated_obs[:check_collision_tv], generated_obs[:check_collision_segment_tv];
		num_iterations=num_rrt_iterations, step_size=step_size,
		radius_goal= .3)
		opt_trajectory_rrt_tv = t-> rrt_nodes_to_path_tv(node_positions, graph_edges, t)[1:end-1]
		opt_trajectory_rrt_tv
	end
	opt_trajectory_rrt_tv = run_rrt_tv(obstacle_setup)
	is_path_valid(opt_trajectory_rrt_tv, obstacle_setup), "smoothness", path_smoothness(opt_trajectory_rrt_tv)
end

# ╔═╡ d24f3e22-0603-11eb-360d-b789b9e18a3a
if n == 2
	q_rrt_tv = plot_obstacles(obstacle_setup, opt_trajectory_rrt_tv)
	PyPlot.plot(eachcol(hcat(opt_trajectory_rrt_tv.(0:.01:1)...)')...)
	q_rrt_tv
end

# ╔═╡ cbaba87a-f876-11ea-0533-756c79f0c51c
begin

results = []
for obs_seed=range_obs_seed
# for solver_seed=range_solver_seed
for speed_obs=range_speed_obs
for n=range_n
	current_run = (obs_seed, solver_seed, n, speed_obs)
	@show (obs_seed, solver_seed, n, speed_obs)
	#@show current_run
	
	obstacle_setup = Dict(:n => n, 
					:a => [-.99 for _ in 1:n],
					:b => [.99 for _ in 1:n],
	:number_obs => number_obs, 
	:speed_obs => speed_obs,
	:radius_obs => radius_obs,
	:obs_seed => obs_seed)

	time_sos = @elapsed run_sos_result =  run_sos(obstacle_setup)
	time_nlp = @elapsed run_nlp_result = run_nlp(obstacle_setup)
	time_rrt = @elapsed run_rrt_result = run_rrt(obstacle_setup)
	time_rrt_tv = @elapsed run_rrt_tv_result = run_rrt_tv(obstacle_setup)

	result_run = Dict(:sos => Dict(:time=>time_sos, :path => run_sos_result),
					:nlp => Dict(:time => time_nlp, :path => run_nlp_result),
					:rrt => Dict(:time => time_rrt, :path => run_rrt_result),
					:rrt_tv => Dict(:time => time_rrt_tv, :path => run_rrt_tv_result),
					)
	push!(results, (obstacle_setup, result_run))		
end
end
end
end


# ╔═╡ 24ecf380-05f3-11eb-2a71-5d9796f7145d
begin
	result_stat = []
	
	println("Evaluating quality of paths...")
	for run=results
		obstacle_setup = run[1]
		@show run[1]
		push!(result_stat, 
				merge(obstacle_setup,
				[Dict(
				Symbol(:valid, "_", method) => is_path_valid( run[2][method][:path], obstacle_setup),
				Symbol(:length, "_", method) => length_path( run[2][method][:path]),
				Symbol(:time, "_", method) =>  run[2][method][:time],
				Symbol(:smooth, "_", method) =>  path_smoothness(run[2][method][:path]),
					)
				for method=(:sos, :rrt, :nlp, :rrt_tv)]...
			))
	end
	println("Evaluating quality of paths done.")
	result_stat = DataFrame(result_stat)
end

# ╔═╡ 081fd8ba-f883-11ea-2dd2-d3d39e7fcf8b
md"# Save to file"

# ╔═╡ 2992287a-0820-11eb-0f1a-bb975f4c110a
result_stat[:, :smooth_rrt_tv]

# ╔═╡ 9f54ed10-f881-11ea-0db8-23fa207bd876
begin
	CSV.write("csv/results_benchmark-$(Dates.now()).csv", result_stat)
	CSV.write("csv/results_benchmark-latest.csv", result_stat)
end

# ╔═╡ Cell order:
# ╟─c3d41f9e-f7de-11ea-2a56-b76588d6ef66
# ╟─79575ad0-f7de-11ea-2e97-97d0955e0194
# ╟─98a3f246-05e8-11eb-201d-a133f7585b7c
# ╠═a39d4a4a-05e7-11eb-2f92-e118170808f8
# ╟─abd44b0a-05e7-11eb-07ed-c7c73096570e
# ╠═b1917388-05e7-11eb-32fb-45c9ace941de
# ╠═ba98a5aa-05e7-11eb-06c9-653431e3ec70
# ╠═c793b196-05e7-11eb-3fd5-dd7fd6b48e49
# ╠═e956f92a-05ea-11eb-3849-c9c3cdf5fb47
# ╟─cb711d12-05e7-11eb-3056-5fa4b9853578
# ╟─cf3ed704-05e7-11eb-042c-bb774ba8306e
# ╟─d47719ac-05e7-11eb-35b4-ed887c5c643f
# ╠═d86e55c8-05e7-11eb-037d-894ea2b56788
# ╠═de847bf6-05e7-11eb-39ca-874f2cdd74fc
# ╟─587fb63e-066d-11eb-3499-0ff6ff490bb2
# ╟─1d3b7e12-05e8-11eb-3f76-3b7dc80b948f
# ╠═e856be64-05e7-11eb-2260-77c2b187e2ed
# ╠═ec1ce6ea-05e7-11eb-3816-11b297344745
# ╠═f06b8128-05e7-11eb-23db-5fbe3f7cd461
# ╠═90161a4a-f7de-11ea-186d-972892ea3c26
# ╟─e590ec80-f87d-11ea-153d-5d9bed4819a5
# ╠═1631c244-05f1-11eb-330a-4b2bdc5073a0
# ╟─2011598c-05ec-11eb-1890-850db5315bd5
# ╠═4d5b9cae-0605-11eb-2cf3-576448d0439e
# ╟─691cc1b2-0604-11eb-01a1-e9db3148a92a
# ╟─c6c0c038-0600-11eb-1c80-b51afcbbab50
# ╟─cb6f6e86-0600-11eb-33ec-416d9a847dc3
# ╟─d358181e-0600-11eb-1ad9-4316519b7784
# ╟─8981d24e-f874-11ea-0685-73fd3c241512
# ╟─bd7b97e4-f7de-11ea-096f-27a885a176c7
# ╠═b3b554bc-05e8-11eb-30ee-13748ef82f3b
# ╠═3eda12e4-0616-11eb-3791-03110d0f8432
# ╟─725babbe-f874-11ea-2552-9d3968d8b8d9
# ╟─72811ea2-f875-11ea-2c97-b1b5aa72c872
# ╟─771b7552-f875-11ea-3b2c-d7bb19dac4f1
# ╟─09836820-05ee-11eb-08df-993638bce61f
# ╟─937f9bac-f874-11ea-1d59-77ff8d0aaef0
# ╠═8bf7a05a-0602-11eb-0940-1b290811b69d
# ╟─28bd2dd0-05ed-11eb-3389-4dc1b4c7bfe1
# ╟─68ed372c-f875-11ea-1e3d-07d9438dc9d5
# ╟─b696d130-05e9-11eb-374c-1b166d28069d
# ╟─2bf3de7c-05ed-11eb-2189-9b7d947b3c19
# ╠═3408f502-05ed-11eb-1316-cdd8235769bf
# ╟─d24f3e22-0603-11eb-360d-b789b9e18a3a
# ╟─43e4f67a-05e4-11eb-29d3-83cd2faaee75
# ╟─6bce6542-05e5-11eb-1322-81821b08c030
# ╟─ebc7c88c-05e6-11eb-13cd-455b20227590
# ╟─fce9d6a4-0603-11eb-10c7-27fc653d8c99
# ╟─ced25ee2-f874-11ea-3db4-bfd63ce72a87
# ╠═5473ac7a-f8ae-11ea-1a13-41a2ebd0b84f
# ╠═cbaba87a-f876-11ea-0533-756c79f0c51c
# ╠═24ecf380-05f3-11eb-2a71-5d9796f7145d
# ╟─8be8ca2e-f87b-11ea-383c-a75157b2cf35
# ╟─69813a14-0615-11eb-2398-31b3f1b9c905
# ╟─081fd8ba-f883-11ea-2dd2-d3d39e7fcf8b
# ╠═2992287a-0820-11eb-0f1a-bb975f4c110a
# ╠═9f54ed10-f881-11ea-0db8-23fa207bd876
