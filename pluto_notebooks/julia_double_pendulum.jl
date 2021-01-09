### A Pluto.jl notebook ###
# v0.12.2

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

# ╔═╡ 3b4cdb9a-0a74-11eb-2144-650c34ba3b05
begin
	using Pkg
	Pkg.activate("..")	
	using Revise
	using Convex, SCS
	using PathPlanningSOS
	using MosekTools
	using JuMP
	using PlutoUI
	using Plots
	using Random
	using PyPlot
	using LinearAlgebra
	using DataFrames
	using LightGraphs
	using KNITRO
	using Statistics
	using MathOptInterface
	using DynamicPolynomials
	md"""A bunch of imports"""
end

# ╔═╡ f282076e-0b63-11eb-2dac-934fcc3518c9
using DelimitedFiles

# ╔═╡ 253d1fb8-0a74-11eb-13a4-5ff45d6ff0d8
md"# Double Pendulum"

# ╔═╡ 34b18f54-0a76-11eb-0999-c7177895da0b
md"t $(@bind plot_at_time Slider(0:.01:1, show_value=true, default=0.))"

# ╔═╡ 4438f574-0a74-11eb-3161-67626572c265
md"""
# Data of the problem

Pos pendulum 1 $(@bind pend_1_x Slider(-1:.01:1, show_value=true, default=-1)) $(@bind pend_1_y Slider(-1:.01:1, show_value=true, default=0.))


Pos pendulum 2 $(@bind pend_2_x Slider(-1:.01:1, show_value=true, default=0.5)) $(@bind pend_2_y Slider(-1:.01:1, show_value=true, default=0.))



number iterations:
$(@bind num_iterations Slider(1:100, show_value=true, default=100))

number pieces:
$(@bind num_pieces Slider(1:100, show_value=true, default=3))


Trade off between minimizing length and rank
$(@bind weight_lenght Slider(0.00:.0001:.001, show_value=true, default=.001))


arm length ``\ell``:
$(@bind ℓ Slider(0.00:.01:1, show_value=true, default=.5))

solver seed 32
$(@bind solver_seed Slider(1:100, show_value=true, default=32))

"""

# ╔═╡ 6db1b048-0ab9-11eb-34f8-3744af737fd7
solver_seed

# ╔═╡ 47496634-0a77-11eb-129d-4b7f5457c280
struct Pendulum
    base::Array{Number,1}
	ℓ::Number
	θ::Array{Number,1}
end

# ╔═╡ 56c80608-0aa9-11eb-0793-eb915a75ce3f
θ_init = [1; 1]

# ╔═╡ 5faa5406-0aa9-11eb-1c6c-1d5b9cb9b1ec
α_init = [π;#3π/4;
	3π/2]

# ╔═╡ 86311590-0aa9-11eb-0cd0-1fc88e36bb7d
θ_final = [0; 0]

# ╔═╡ 91d89bf4-0aa9-11eb-2aa1-cbef22468866
α_final = [α_init[1], 3π/2-2]

# ╔═╡ b794aef8-0a79-11eb-3f2f-77bb5bebad51
md"""
``\theta`` $(@bind θ_1 Slider(0:.01:π/2, show_value=true, default=θ_init[1] )) $(@bind θ_2 Slider(-π/2:.01:π/2, show_value=true, default=θ_init[2]))

``\alpha`` $(@bind α_1 Slider(3π/4:π:π, show_value=true, default=α_init[1])) $(@bind α_2 Slider(π/2:.01:3π/2, show_value=true, default=α_init[2]))

"""

# ╔═╡ 60101a3e-0a77-11eb-019a-f3c7b7ace024
begin
	pend_1 = Pendulum([pend_1_x, pend_1_y], ℓ, [θ_1, θ_2])
	pend_2 = Pendulum([pend_2_x, pend_2_y], ℓ, [α_1, α_2])
	pend_1, pend_2
end

# ╔═╡ e5f31194-0a7a-11eb-153b-430032c4340c
function get_pend_segments(pend)
	middle_joint = pend.base .+ ℓ .* [cos(pend.θ[1]), sin(pend.θ[1])]
	end_joint = middle_joint .+ ℓ.*  [cos(pend.θ[2]), sin(pend.θ[2])]
	joints = [pend.base,
				middle_joint,
				end_joint]
end

# ╔═╡ a1cbfa84-0a76-11eb-29e4-7d1f1573d1c2
function plot_pend(pend; c="red", plot_ref_line=true)
	joints = get_pend_segments(pend)
	
	P = hcat(joints...)
	PyPlot.plot(P[1, :], P[2, :], c=c)
	PyPlot.scatter(P[1, :], P[2, :], s=10, c=c)
	
	if plot_ref_line
	for joint ∈ joints[1:2]
		PyPlot.plot([joint[1], joint[1]+ℓ], [joint[2], joint[2]], c="k", ls="--")
	end
	end
end

# ╔═╡ 1464189e-0a7a-11eb-2a89-9148d08af914
md"# Collision checker"

# ╔═╡ 1b32b464-0a7a-11eb-2e22-27d4024f85ef
function segments_intersect(s_1, s_2)
	"""
	Code from here: https://stackoverflow.com/a/9997374/440221
	"""
	ccw(A,B,C) = (C[2]-A[2]) * (B[1]-A[1]) > (B[2]-A[2]) * (C[1]-A[1])
	A, B = s_1
	C, D = s_2
	(ccw(A,C,D) != ccw(B,C,D)) & (ccw(A,B,C) != ccw(A,B,D))
end

# ╔═╡ d6234202-0a7a-11eb-0578-75877aea6517
function pend_intersect(pend_1, pend_2)
	joints_1 = get_pend_segments(pend_1)
	joints_2 = get_pend_segments(pend_2)
	
	for s_1 in zip(joints_1[1:end-1], joints_1[2:end])
		for s_2 in zip(joints_2[1:end-1], joints_2[2:end])
			if segments_intersect(s_1, s_2)
				return true
			end
		end
	end
	false
end

# ╔═╡ fb78d7be-0a7b-11eb-111b-57cd4801ae28
md"# Obstacle Map computation"

# ╔═╡ 2ef39aac-0a7c-11eb-374e-6f481389bb68
md"""

``N``  $(@bind N Slider(1:50, show_value=true, default=40))
"""

# ╔═╡ 45d86514-0aa6-11eb-1aec-f777eaec430b
ranges = [
		(0 -π/2, π/2 +π/2, N),
		(-π/2 -π/2, π/2+π/2, N),
		(3pi/4, 3pi/4, 1),
		#(1. *π, 1. *π, 1),
		(π -π/2, 3π/2+π/2, N),
		]

# ╔═╡ 2d6f5bc6-0a7c-11eb-198d-67716bac325e
begin
	
	angle_mesh = collect(Iterators.product([range(θmin, stop=θmax, length=N) for (θmin,θmax,N)=ranges]...))
	md"Angle mesh size: $(prod(size(angle_mesh)))"
end

# ╔═╡ c289057c-0a7c-11eb-1dfa-17300bb73af8
begin
	collision_map = []
	for (θ_1, θ_2, α_1, α_2) in angle_mesh
		pend_1 = Pendulum([pend_1_x, pend_1_y], ℓ, [θ_1, θ_2])
		pend_2 = Pendulum([pend_2_x, pend_2_y], ℓ, [α_1, α_2])
		push!(collision_map, (θ_1, θ_2, α_1, α_2, pend_intersect(pend_1, pend_2)))
	end
	collision_map = DataFrame(collision_map)
	names!(collision_map, Symbol.(["θ_1", "θ_2", "α_1", "α_2", "collide"]))
	describe(collision_map)
end

# ╔═╡ 2420b20a-0a80-11eb-322a-758db1d10e76
begin
	configs_that_collide = filter(row -> row[:collide], collision_map)
	describe(configs_that_collide)
end

# ╔═╡ dd1d1afe-0a81-11eb-08c2-2fc13dcf722c
begin
	Plots.plot(configs_that_collide.θ_1,
		configs_that_collide.θ_2,
		configs_that_collide.α_2,
		seriestype=:scatter, markersize = 7
		)
	Plots.xlims!( ranges[1][1:end-1]...)
	Plots.ylims!( ranges[2][1:end-1]...)
	Plots.zlims!( ranges[4][1:end-1]...)
	
	
	Plots.xlims!( -π, π)
	Plots.ylims!( -π, π)
	Plots.zlims!(0, 2π)
end

# ╔═╡ 112d723c-0a8d-11eb-0cec-ab7a9c30a3bd
function fit_sphere(points)
	dist(a, b) = LinearAlgebra.norm(a .- b)
	center = Statistics.mean(points, dims=1)
	radius = maximum([dist(center, points[i, :]) for i=1:size(points,1)])
	
	center, radius
end

# ╔═╡ fee139e8-0a8d-11eb-23c4-f1fdaf19b811
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# ╔═╡ d3248326-0a8d-11eb-270c-99e09c28bdcb
function fit_ellipse(points; ε=1e-2)
	# { x | x = C*u + d, ||u||_2 <= 1 }
	
	
	n = size(points, 2)
	if size(points, 1) == 0
		return zeros(n, n), vec(zeros(n, 1))
	end
	
	
	P = Variable(n, n)
	q = Variable(n)
	problem = maximize(logdet(P))
	
	for i=1:size(points, 1)
		xi = points[i, :]
		problem.constraints += norm(P*xi - q) <= 1 - ε
	end
	Convex.solve!(problem, SCS.Optimizer)
	P, q = P.value, q.value
	C = inv(P)
	d = vec(C * q) 
	C, d
end

# ╔═╡ 24fd89fc-0a8d-11eb-0628-018cc5b2df9a
points_to_fit = hcat(configs_that_collide.θ_1,
		configs_that_collide.θ_2,
		configs_that_collide.α_2)

# ╔═╡ ee8e79d8-0a8e-11eb-3bc3-6b8e8d61e74c
function plot_ellipse(C::Array{Float64, 2}, d::Array{Float64, 1}; ngrid=20)
	"""
	Plots the set { C*u + d | u in R^n }
	"""
	# Set of all spherical angles:
	u = range(0,stop=2*π,length=ngrid);
	v = range(0,stop=π,length=ngrid);

	x = cos.(u) * sin.(v)';
	y = sin.(u) * sin.(v)';
	z = ones(ngrid) * cos.(v)';
	X, Y, Z = vec(x), vec(y), vec(z)
	P = hcat(X, Y, Z) * C
	X, Y, Z = P[:, 1], P[:, 2], P[:, 3]
	
    X, Y, Z = X .+ d[1], Y .+ d[2], Z .+ d[3]
	scatter!(X, Y, Z, aspect_ratio=:equal, label="", markersize=2)
	
end

# ╔═╡ e0821706-0aa6-11eb-0c7b-53986e8479ec
function plot_cube(range_x, range_y, range_z)
	for i=1:2
		Plots.plot!([range_x[i],range_x[i],range_x[i],range_x[i],range_x[i]],
					[range_y[1],range_y[2],range_y[2],range_y[1],range_y[1]],
					[range_z[1],range_z[1],range_z[2],range_z[2],range_z[1]],
					c="black", label="", lw=4
					)
		for j=1:2
		Plots.plot!([range_x[1],range_x[2]],
		[range_y[i],range_y[i]],
		[range_z[j],range_z[j]], c="black", label="", lw=4
		)
		end
	end

end

# ╔═╡ 9786dfe4-0a94-11eb-03f9-459ad459e1b5
md"Camera $(@bind cam_angle Slider(0:1:100, show_value=true, default=20.))

Camera $(@bind cam_angle2 Slider(0:1:100, show_value=true, default=40.))

"

# ╔═╡ ae951cb8-0a90-11eb-302f-9f4996453a8f
begin
	C, d = fit_ellipse(points_to_fit)
	Plots.plot(configs_that_collide.θ_1,
		configs_that_collide.θ_2,
		configs_that_collide.α_2,
		seriestype=:scatter, markersize = 2,
		label="",
		camera = (cam_angle, cam_angle2)
		)
	#Plots.plot(camera = (cam_angle, cam_angle2))
	plt = plot_ellipse(C, d)

	X_init_final = hcat([θ_init..., α_init[2]], [θ_final..., α_final[2]],)
	xinit = X_init_final[:, 1]
	xfinal = X_init_final[:, 2]
	Plots.scatter!([X_init_final[i,:] for i=1:3]..., label="x init final")
	Plots.xlims!( ranges[1][1:end-1]...)
	Plots.ylims!( ranges[2][1:end-1]...)
	Plots.zlims!( ranges[4][1:end-1]...)
	
	#plot_cube(ranges[1], ranges[2], ranges[4])

	
# 	Plots.xlims!( -π, π)
# 	Plots.ylims!( -π, π)
# 	Plots.zlims!(0, 2π)
	

end

# ╔═╡ f41d6656-0a7b-11eb-04b1-47b6d68ba35e
md"# Motion planning"

# ╔═╡ cda005a6-0ab5-11eb-335d-f9f8a9f351e3
ε = 1e-3

# ╔═╡ 7461d972-0aaa-11eb-092d-110ff7b2a756
P = inv(C + ε*I) * (inv(C+ ε*I))'

# ╔═╡ 1c6cf804-0aa6-11eb-390e-17f84902f79e
obstacles = [
 	(t, x) -> (x .- d)' * P * (x .- d) - 1,
	[(t, x) -> x[i] - ri[1] for (i,ri)=zip(1:3, (ranges[1], ranges[2], ranges[4]))]...,
	[(t, x) -> ri[2] - x[i] for (i,ri)=zip(1:3, (ranges[1], ranges[2], ranges[4]))]...,
	]

# ╔═╡ f78341ee-0a7b-11eb-3afe-31e33b156a76
begin
	# compute optimal piece-wise linear trajectory
	Random.seed!(solver_seed)
	n = 3
	deg = 2
	world_x = 2π
	opt_trajectory = find_path_using_heuristic(n, obstacles, world_x, 
		xinit, xfinal,
	    deg, num_pieces, solver,
	    weight_lenght,
	    num_iterations,
	    seed=solver_seed)
	md"Computing path"
end

# ╔═╡ a96bfb76-0aae-11eb-1851-0dba02fb43e9
begin
	intersection_times = []
	for t=0:.001:1
		x_t = opt_trajectory(t)
		pend_1_t = Pendulum([pend_1_x, pend_1_y], ℓ, x_t[1:2])
		pend_2_t = Pendulum([pend_2_x, pend_2_y], ℓ, [α_1, x_t[3]])
		if pend_intersect(pend_1_t, pend_2_t)
			push!(intersection_times, t)
		end
	end
	intersection_times
end

# ╔═╡ 016bf988-0aab-11eb-3694-d304b5776068
[f(0, xfinal) for f=obstacles]

# ╔═╡ c32f8b52-0ab5-11eb-232e-8b4b5a0345df
P

# ╔═╡ 18c74878-0aae-11eb-3f6a-3395e7f8ce2b
[(t, f(t, opt_trajectory(t)) >= 0) for t=0.5:.01:.6
	for f=obstacles[1:1]]

# ╔═╡ 21e2a0e4-0aae-11eb-2ff1-61643ca2b7e3
 obstacles[1](.53, opt_trajectory(.61)) 

# ╔═╡ 9ea9132e-0ab6-11eb-27ba-ffed48a1b08a
opt_trajectory(.61)

# ╔═╡ 790af880-0ab6-11eb-2aa4-a362fc148651
all([f(t, opt_trajectory(t)) >= 0 for t=0:.01:1 for f=obstacles])

# ╔═╡ 670959f4-0aae-11eb-3b32-cf789d2b0a50
 LinearAlgebra.norm(inv(C) * (opt_trajectory(.53) .- d))

# ╔═╡ 72aa2a68-0aae-11eb-2d31-dd4fbacd751b
(opt_trajectory(.52) .- d)' * P * (opt_trajectory(.53) .- d) - 1

# ╔═╡ 730f3f8a-0abc-11eb-32b7-d78a7c94289f
md"# Plots for paper"

# ╔═╡ 07438a9e-0abd-11eb-3ecc-33556729a66f
function plot_everything_at_time(t, opt_trajectory)
	q = PyPlot.figure()
	x_t = opt_trajectory(t)
	pend_1_t = Pendulum([pend_1_x, pend_1_y], ℓ, x_t[1:2])
	pend_2_t = Pendulum([pend_2_x, pend_2_y], ℓ, [α_1, x_t[3]])
	plot_pend(pend_1_t, c="blue",plot_ref_line=false)
	plot_pend(pend_2_t, c="green",plot_ref_line=false)
	PyPlot.ylim(-.5, 1)
	PyPlot.xlim(-1.2, 1)
	#PyPlot.title("Initial configuration")
	PyPlot.axes().set_aspect("equal")
	PyPlot.tight_layout()
	gca()[:axis]("off")
	q
end

# ╔═╡ 1a2b598c-0aab-11eb-04c3-dd30094ced04
begin
	
	@info "PLotting"
	q_sos = plot_everything_at_time(plot_at_time, opt_trajectory)
	PyPlot.title("t = $plot_at_time | SOS")
	q_sos
end

# ╔═╡ 300b2202-0abd-11eb-04d5-0753ba90a2ce
plot_everything_at_time(0, opt_trajectory)

# ╔═╡ 7e0c70a6-0abc-11eb-3d28-2fcd5b0fdb4d
begin
	qq = plot_everything_at_time(0, opt_trajectory)
	PyPlot.title("Initial configuration")
end

# ╔═╡ 89a69202-0abc-11eb-2d67-47f5fbf662d3
begin
	qT = plot_everything_at_time(1., opt_trajectory)
	PyPlot.title("Goal configuration", fontsize=30)
	qT
end

# ╔═╡ 04a4e1de-0abd-11eb-38ad-51098e49e416
begin
	for t=0:.01:1
		qt = plot_everything_at_time(t, opt_trajectory)
		t_str = string(Int(round(t*100, digits=3)), pad=3)
		PyPlot.title("MMP",  fontsize=30)
		PyPlot.savefig("Imgs/double_pendulum_$(t_str)_sos.png")
	end
end

# ╔═╡ a92afcb6-0d80-11eb-0dbb-514cce26ff82


# ╔═╡ 91e45136-0b73-11eb-229e-8fbeda95fb71
md"# RRT"

# ╔═╡ f18c22f8-0b73-11eb-09b4-f5d16e6baf59
function rrt_find_path(dim, start, goal, check_collision, check_collision_segment, world_ranges; num_iterations=10, step_size=.2, radius_goal=.1)
	
	node_positions = [ start ]
	graph_edges = []
	#Random.seed!(obs_seed)
	
	for counter=1:num_iterations
		
		# pick a random point
		new_x = [ Random.rand(Float64) * (world_ranges[i][2] - world_ranges[i][1]) + world_ranges[i][1] for i=1:dim ]
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


# ╔═╡ f47ba150-0b73-11eb-17cf-a52d804d0d63
function rrt_nodes_to_path(node_positions, graph_edges, t)
	g = LightGraphs.SimpleGraph(size(node_positions, 1))
	for e=graph_edges
		LightGraphs.add_edge!(g, e...)
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


# ╔═╡ 962c2c0a-0b73-11eb-01c4-43167ccff7fe
begin
	function check_collision(q)
		all([f(t, q) >= 0 for t=0:.01:1 for f=obstacles])
	end
	
	function check_collision_segment(q1, q2)
		all([check_collision(α .* q1 .+ (1-α) .* q2) for α=0:.01:1])
	end
	
end

# ╔═╡ 10745848-0b74-11eb-1df9-250d5dc71511
begin
	
	# compute optimal piece-wise linear trajectory
	Random.seed!(solver_seed)
	num_iterations_rrt = 1000
	step_size = 1.
	
	rrt_path_valid, node_positions, graph_edges = rrt_find_path(n,xinit, xfinal, check_collision, check_collision_segment, [ranges[1], ranges[2], ranges[4]]; num_iterations=num_iterations_rrt, step_size=step_size, radius_goal=1.)
	opt_trajectory_rrt = t-> rrt_nodes_to_path(node_positions, graph_edges, t)
	md"Computing path RRT"
end

# ╔═╡ 87e64760-0b8d-11eb-31e3-2db3037716db
begin
	
	@info "PLotting"
	q_rrt = plot_everything_at_time(plot_at_time, opt_trajectory_rrt)
	PyPlot.title("t = $plot_at_time | RRT")
	q_rrt
end

# ╔═╡ aa39fa40-0b8b-11eb-0791-2fd43fe1d849
begin
	for t=0:.01:1
		qt = plot_everything_at_time(t, opt_trajectory_rrt)
		t_str = string(Int(round(t*100, digits=3)), pad=3)
		PyPlot.title("RRT", fontsize=30)
		PyPlot.savefig("Imgs/double_pendulum_$(t_str)_rrt.png")
	end
end

# ╔═╡ 7edcceb8-0b75-11eb-0af0-09a186301228
rrt_path_valid

# ╔═╡ 1a8eb1aa-0b76-11eb-1bc9-37233b758d43
opt_trajectory_rrt(.5)

# ╔═╡ b995a532-0b8b-11eb-0b64-a9bd35b001ce
md"# NLP"

# ╔═╡ c25ab5e0-0b8b-11eb-37f8-9544d40c60d5
function shortest_path_nlp(
		n,
		a, b,
		num_pieces,
		obstacles,
		)
	@info "Model"
	
	model = Model(optimizer_with_attributes(KNITRO.Optimizer, "honorbnds" => 1,
											"outlev" => 0, "algorithm" => 4,
											))
	@variable(model, u[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	@variable(model, v[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	
		lengths = @NLexpression model [i=1:size(v,1)] sqrt(sum(v[i, j]^2 for j=1:n))
		
	
	
		@info "Constraints"

		for i=1:num_pieces
			for t=0:.1:1
				# xt = u[i,:] .+ t .* v[i,:] .- d

				# Ellipse
				@NLconstraint model  sum((u[i,ii] + t * v[i,ii] - d[ii]) * P[ii, jj] * (u[i,jj] + t * v[i,jj] - d[jj])  for ii=1:n for jj=1:n)>= 1.
			end
		end
	
		down, up = [[ranges[1][i], ranges[2][i], ranges[4][i]] for i=1:2]
		for i=1:num_pieces-1
			@constraint model  u[i,:] .>=  down
			@constraint model  u[i,:] .<=  up
		end
		for i=1:num_pieces-1
			@constraint model  u[i,:] .+ v[i,:] .== u[i+1,:]
		end
		
		@constraint model u[1,:] .== a
		@constraint model u[end,:] .+ v[end,:] .== b;
		
		
	
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

# ╔═╡ d4ea39ec-0b8b-11eb-06aa-2d61fc0cbb1f
begin
	
	# compute optimal piece-wise linear trajectory
	Random.seed!(solver_seed)

	_, _, opt_trajectory_nlp = shortest_path_nlp(
		n, 
		xinit, xfinal,
		num_pieces,
		obstacles,
		)
	md"Computing path NLP"
end

# ╔═╡ bc12162c-0b8d-11eb-1210-0736964d4a3b
begin
	
	@info "PLotting"
	q_nlp = plot_everything_at_time(plot_at_time, opt_trajectory_nlp)
	PyPlot.title("t = $plot_at_time | NLP")
	q_nlp
end

# ╔═╡ eafa3134-0aa6-11eb-2b38-5de4cbe8f1c7
begin
	
	Plots.plot(configs_that_collide.θ_1,
		configs_that_collide.θ_2,
		configs_that_collide.α_2,
		seriestype=:scatter, markersize = 1,
		label="",
		camera = (cam_angle, cam_angle2)
		)
	
	plot_cube(ranges[1], ranges[2], ranges[4])
	plot_ellipse(C, d)
	Plots.scatter!( [ [xinit[i]] for i=1:3]..., label="xinit")
	Plots.scatter!( [ [xfinal[i]] for i=1:3]..., label="xfinal")
	opt_path_discretized = hcat(opt_trajectory.(0:.01:1)...)
	path_valid = all([f(t, opt_trajectory(t)) >= 0 for t=0:.01:1 for f=obstacles])
	Plots.plot!([opt_path_discretized[i, :] for i=1:size(opt_path_discretized,1)]...,
	lw=4, c="blue")
	
	# RRT
	opt_path_rrt_discretized = hcat(opt_trajectory_rrt.(0:.01:1)...)
	Plots.plot!([opt_path_rrt_discretized[i, :] for i=1:size(opt_path_rrt_discretized,1)]...,
	lw=4, c="red")
	
		# RRT
	opt_path_nlp_discretized = hcat(opt_trajectory_nlp.(0:.01:1)...)
	Plots.plot!([opt_path_nlp_discretized[i, :] for i=1:size(opt_path_nlp_discretized,1)]...,
	lw=4, c="green")
	
	
	Plots.xlabel!("θ_1")
	Plots.ylabel!("θ_2")
	#Plots.zlabel!("α_1")
	#Plots.title!("Path valid? $(path_valid)")
	Plots.title!("Obstacle in joint space")
end

# ╔═╡ 5ecd3cc0-0b8e-11eb-1c00-f563b7c2e53e
begin
	for t=0:.01:1
		qt = plot_everything_at_time(t, opt_trajectory_nlp)
		t_str = string(Int(round(t*100, digits=3)), pad=3)
		PyPlot.title("NLP", fontsize=30)
		PyPlot.savefig("Imgs/double_pendulum_$(t_str)_nlp.png")
	end
end

# ╔═╡ 589df188-0b8d-11eb-21d5-19354c7eeea0
opt_trajectory_nlp(.5)

# ╔═╡ d7d92636-0b63-11eb-2b57-e1a2305be424
md"# Export"

# ╔═╡ db94cf78-0b63-11eb-0788-6d52f480f75c
writedlm( "csv/obstacles_in_joint_space.csv",  points_to_fit, ',')

# ╔═╡ fed0ae3e-0b69-11eb-061c-afd56375467a
begin
	function export_ellipse2(C::Array{Float64, 2}, d::Array{Float64, 1}; ngrid=20)
		"""
		Plots the set { C*u + d | u in R^n }
		"""
		# Set of all spherical angles:
		u = range(0,stop=2*π,length=ngrid);
		v = range(0,stop=π,length=ngrid);
	
		x = cos.(u) * sin.(v)';
		y = sin.(u) * sin.(v)';
		z = ones(ngrid) * cos.(v)';
		X, Y, Z = vec(x), vec(y), vec(z)
		P = hcat(X, Y, Z) * C
		X, Y, Z = P[:, 1], P[:, 2], P[:, 3]
		
	    X, Y, Z = X .+ d[1], Y .+ d[2], Z .+ d[3]
		hcat(vec(X), vec(Y), vec(Z))
	end
	writedlm( "csv/fittedellipse.csv",  export_ellipse2(C, d; ngrid=50), ",")
end

# ╔═╡ 6d358474-0b7b-11eb-3548-9b287d716ea6
opt_trajectory_nlp(1)

# ╔═╡ 4f45e542-0b78-11eb-3e98-458ace58062b
begin
	writedlm( "csv/optsostrajectory.csv",  hcat(opt_trajectory.(0:.05:1)...)', ',')
	writedlm( "csv/optrrttrajectory.csv",  hcat(opt_trajectory_rrt.(0:.05:1)...)', ',')
	writedlm( "csv/optnlptrajectory.csv",  hcat(opt_trajectory_nlp.(0:.05:1)...)', ',')
	
end

# ╔═╡ ef3aaa1a-0d80-11eb-0877-3b0404ffe0b0
writedlm( "csv/opttrajectoryall.csv", vcat([ hcat(traj.(0:.01:1)...) for traj in (opt_trajectory, opt_trajectory_rrt, opt_trajectory_nlp)]...)', ',')

# ╔═╡ 03a4080e-0bed-11eb-24b5-f320fe844c26
md"""
```cp  csv/opt*csv csv/fittedellipse.csv csv/obstacles_in_joint_space.csv ~/Overleaf/Path_Planning_ICRA_2020/Imgs/```

```
cp Imgs/*png  ~/Overleaf/Path_Planning_ICRA_2020/Imgs/double_pendulum```
"""

# ╔═╡ 3b25ae12-0bed-11eb-2ebc-ff7c95c307f5
get_pend_segments(Pendulum([pend_1_x, pend_1_y], ℓ, xinit[1:2]))

# ╔═╡ 02567d96-0bf5-11eb-0e34-5941a842ba1e
get_pend_segments(Pendulum([pend_2_x, pend_2_y], ℓ, [α_1, xinit[3]]))

# ╔═╡ c0ab9eae-0bf8-11eb-3fd4-af72b9fa68e9
xinit

# ╔═╡ Cell order:
# ╠═253d1fb8-0a74-11eb-13a4-5ff45d6ff0d8
# ╠═3b4cdb9a-0a74-11eb-2144-650c34ba3b05
# ╠═6db1b048-0ab9-11eb-34f8-3744af737fd7
# ╟─34b18f54-0a76-11eb-0999-c7177895da0b
# ╟─1a2b598c-0aab-11eb-04c3-dd30094ced04
# ╟─87e64760-0b8d-11eb-31e3-2db3037716db
# ╟─bc12162c-0b8d-11eb-1210-0736964d4a3b
# ╠═eafa3134-0aa6-11eb-2b38-5de4cbe8f1c7
# ╠═a96bfb76-0aae-11eb-1851-0dba02fb43e9
# ╠═4438f574-0a74-11eb-3161-67626572c265
# ╠═47496634-0a77-11eb-129d-4b7f5457c280
# ╠═56c80608-0aa9-11eb-0793-eb915a75ce3f
# ╠═5faa5406-0aa9-11eb-1c6c-1d5b9cb9b1ec
# ╠═86311590-0aa9-11eb-0cd0-1fc88e36bb7d
# ╠═91d89bf4-0aa9-11eb-2aa1-cbef22468866
# ╟─b794aef8-0a79-11eb-3f2f-77bb5bebad51
# ╠═60101a3e-0a77-11eb-019a-f3c7b7ace024
# ╟─e5f31194-0a7a-11eb-153b-430032c4340c
# ╠═a1cbfa84-0a76-11eb-29e4-7d1f1573d1c2
# ╟─1464189e-0a7a-11eb-2a89-9148d08af914
# ╠═1b32b464-0a7a-11eb-2e22-27d4024f85ef
# ╟─d6234202-0a7a-11eb-0578-75877aea6517
# ╟─fb78d7be-0a7b-11eb-111b-57cd4801ae28
# ╟─2ef39aac-0a7c-11eb-374e-6f481389bb68
# ╠═45d86514-0aa6-11eb-1aec-f777eaec430b
# ╟─2d6f5bc6-0a7c-11eb-198d-67716bac325e
# ╟─c289057c-0a7c-11eb-1dfa-17300bb73af8
# ╟─2420b20a-0a80-11eb-322a-758db1d10e76
# ╟─dd1d1afe-0a81-11eb-08c2-2fc13dcf722c
# ╟─112d723c-0a8d-11eb-0cec-ab7a9c30a3bd
# ╠═fee139e8-0a8d-11eb-23c4-f1fdaf19b811
# ╠═d3248326-0a8d-11eb-270c-99e09c28bdcb
# ╠═24fd89fc-0a8d-11eb-0628-018cc5b2df9a
# ╠═ee8e79d8-0a8e-11eb-3bc3-6b8e8d61e74c
# ╟─e0821706-0aa6-11eb-0c7b-53986e8479ec
# ╟─9786dfe4-0a94-11eb-03f9-459ad459e1b5
# ╠═ae951cb8-0a90-11eb-302f-9f4996453a8f
# ╟─f41d6656-0a7b-11eb-04b1-47b6d68ba35e
# ╠═cda005a6-0ab5-11eb-335d-f9f8a9f351e3
# ╠═7461d972-0aaa-11eb-092d-110ff7b2a756
# ╠═1c6cf804-0aa6-11eb-390e-17f84902f79e
# ╠═f78341ee-0a7b-11eb-3afe-31e33b156a76
# ╠═016bf988-0aab-11eb-3694-d304b5776068
# ╟─c32f8b52-0ab5-11eb-232e-8b4b5a0345df
# ╠═18c74878-0aae-11eb-3f6a-3395e7f8ce2b
# ╠═21e2a0e4-0aae-11eb-2ff1-61643ca2b7e3
# ╠═9ea9132e-0ab6-11eb-27ba-ffed48a1b08a
# ╠═790af880-0ab6-11eb-2aa4-a362fc148651
# ╠═670959f4-0aae-11eb-3b32-cf789d2b0a50
# ╠═72aa2a68-0aae-11eb-2d31-dd4fbacd751b
# ╠═730f3f8a-0abc-11eb-32b7-d78a7c94289f
# ╠═07438a9e-0abd-11eb-3ecc-33556729a66f
# ╠═300b2202-0abd-11eb-04d5-0753ba90a2ce
# ╠═7e0c70a6-0abc-11eb-3d28-2fcd5b0fdb4d
# ╠═89a69202-0abc-11eb-2d67-47f5fbf662d3
# ╠═04a4e1de-0abd-11eb-38ad-51098e49e416
# ╠═aa39fa40-0b8b-11eb-0791-2fd43fe1d849
# ╠═5ecd3cc0-0b8e-11eb-1c00-f563b7c2e53e
# ╠═a92afcb6-0d80-11eb-0dbb-514cce26ff82
# ╟─91e45136-0b73-11eb-229e-8fbeda95fb71
# ╟─f18c22f8-0b73-11eb-09b4-f5d16e6baf59
# ╟─f47ba150-0b73-11eb-17cf-a52d804d0d63
# ╟─962c2c0a-0b73-11eb-01c4-43167ccff7fe
# ╟─10745848-0b74-11eb-1df9-250d5dc71511
# ╠═7edcceb8-0b75-11eb-0af0-09a186301228
# ╠═1a8eb1aa-0b76-11eb-1bc9-37233b758d43
# ╟─b995a532-0b8b-11eb-0b64-a9bd35b001ce
# ╟─c25ab5e0-0b8b-11eb-37f8-9544d40c60d5
# ╠═d4ea39ec-0b8b-11eb-06aa-2d61fc0cbb1f
# ╠═589df188-0b8d-11eb-21d5-19354c7eeea0
# ╠═d7d92636-0b63-11eb-2b57-e1a2305be424
# ╠═f282076e-0b63-11eb-2dac-934fcc3518c9
# ╠═db94cf78-0b63-11eb-0788-6d52f480f75c
# ╠═fed0ae3e-0b69-11eb-061c-afd56375467a
# ╠═6d358474-0b7b-11eb-3548-9b287d716ea6
# ╠═4f45e542-0b78-11eb-3e98-458ace58062b
# ╠═ef3aaa1a-0d80-11eb-0877-3b0404ffe0b0
# ╠═03a4080e-0bed-11eb-24b5-f320fe844c26
# ╠═3b25ae12-0bed-11eb-2ebc-ff7c95c307f5
# ╠═02567d96-0bf5-11eb-0e34-5941a842ba1e
# ╠═c0ab9eae-0bf8-11eb-3fd4-af72b9fa68e9
