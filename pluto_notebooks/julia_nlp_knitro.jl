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

# ╔═╡ f1cc2fae-05df-11eb-1656-8d20206e0579
begin
	using Pkg
	Pkg.activate("..")
	using JuMP
	using PlutoUI
	using Plots
	using Random
	using PyPlot
	using LinearAlgebra
	# using Ipopt
	using DynamicPolynomials
	using KNITRO
	using PathPlanningSOS
	# using MultivariateMoments
	md"""A bunch of imports"""
end

# ╔═╡ 066e1fda-05e0-11eb-08e4-938847b33bb0
dist_squared(a, b) = sum((a .- b).^2)

# ╔═╡ 14751a10-05e1-11eb-2349-5d82d455daae
n = 2

# ╔═╡ 2952fd22-05e0-11eb-0582-c7c4f79c294b

md"""
# Data of the problem

world size:

$(@bind world_x Slider(0.00:.01:1, show_value=true, default=1.)) $(@bind world_y Slider(0.00:.01:1, show_value=true, default=1.))

start: 

$(@bind start_x Slider(-1:.01:1, show_value=true, default=-.99)) $(@bind start_y Slider(-1:.01:1, show_value=true, default=-.98))

goal:

$(@bind end_x Slider(0.00:.01:1, show_value=true, default=1.)) $(@bind end_y Slider(0.00:.01:1, show_value=true, default=1.))


Number obstacles:
$(@bind number_obs Slider(1:100, show_value=true, default=10))


Number pieces:
$(@bind num_pieces Slider(1:10, show_value=true, default=3))



Radius of obstacles:
$(@bind radius_obs Slider(0.00:.01:1, show_value=true, default=.1))



Speed obstacles
$(@bind speed_obs Slider(0.00:.1:1, show_value=true, default=0.))


obstacle seed:
$(@bind obs_seed Slider(1:100, show_value=true, default=0))

"""

# ╔═╡ 2da65c84-05e0-11eb-3d50-f78137db8304
function shortest_path_nlp(a, 
			b,
		num_pieces,
		edge_size,
		obs_pos, obs_vel, radius_obs)
	@info "Model"
	
	model = Model(optimizer_with_attributes(KNITRO.Optimizer, "honorbnds" => 1,
											"outlev" => 1, "algorithm" => 4))
	@variable(model, u[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	@variable(model, v[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	
	@info "Constraints"
		lengths = @NLexpression model [i=1:size(v,1)] sqrt(sum(v[i, j]^2 for j=1:n))
		
		for j=size(obs_pos, 1)
			for i=1:num_pieces
				for t=0:.1:1
					obs_pos_t = obs_pos[j, :] .+ t .* speed_obs .* obs_vel[j, :]
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

# ╔═╡ 5e1e7868-05e0-11eb-2c94-e7cf4110bfe0
begin
	Random.seed!(obs_seed)
	obs_pos = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obs_vel = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obstacles = [
    	(t, x) -> sum( (x .- (obs_pos[i, :] .+ t .* speed_obs .* obs_vel[i, :])).^2 ) - radius_obs^2 for i=1:number_obs
	]
		
	md"obstacles ($number_obs)  generated of radius $(radius_obs)"
end

# ╔═╡ c43eb4c6-05e1-11eb-18aa-71dd85ef61a1
 opt_value, opt_trajectory_length, opt_trajectory = shortest_path_nlp([start_x, start_y], 
	[end_x, end_y],
		num_pieces,
	world_x,
		obs_pos, obs_vel, radius_obs)

# ╔═╡ 8733fa58-05e1-11eb-2ffb-5389caf1b221
md"""# Plot

t $(@bind plot_at_time Slider(0:.01:1, show_value=true, default=0.))"""

# ╔═╡ 751d5b0c-05e1-11eb-2a59-314e6476747f
begin
q = PyPlot.figure()
	PathPlanningSOS.plot_at_time(plot_at_time, world_x, [start_x, start_y], [end_x, end_y],
						obstacles, opt_trajectory)
	PyPlot.title("t = $plot_at_time")
	PyPlot.axes().set_aspect("equal")
	q
end

# ╔═╡ Cell order:
# ╟─f1cc2fae-05df-11eb-1656-8d20206e0579
# ╟─066e1fda-05e0-11eb-08e4-938847b33bb0
# ╠═14751a10-05e1-11eb-2349-5d82d455daae
# ╟─2952fd22-05e0-11eb-0582-c7c4f79c294b
# ╟─2da65c84-05e0-11eb-3d50-f78137db8304
# ╟─c43eb4c6-05e1-11eb-18aa-71dd85ef61a1
# ╟─5e1e7868-05e0-11eb-2c94-e7cf4110bfe0
# ╟─8733fa58-05e1-11eb-2ffb-5389caf1b221
# ╠═751d5b0c-05e1-11eb-2a59-314e6476747f
