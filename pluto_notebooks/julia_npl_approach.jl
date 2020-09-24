### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 89490f12-f8b3-11ea-267f-df907a3b7725
begin
	using Pkg
	Pkg.activate("..")
	Pkg.add("Ipopt")
	#Pkg.add(PackageSpec(url="https://github.com/bachirelkhadir/PathPlanningSOS.jl"))
	using PathPlanningSOS
	using MosekTools
	using JuMP
	using PlutoUI
	using Plots
	using Random
	using PyPlot
	using LinearAlgebra
	using Ipopt
	using DynamicPolynomials
	using MultivariateMoments
	md"""A bunch of imports"""
end

# ╔═╡ 06250748-f8b4-11ea-2dbc-492b649a4fa2

## Helper functions
dist_squared(a, b) = sum((a .- b).^2)



# ╔═╡ 08c04ee0-f8b4-11ea-0ae4-ebe7288fe0bf
function plot_levelset(g, x, y)
	min_x, max_x = PyPlot.xlim()
	min_y, max_y = PyPlot.ylim()
	xs = collect(min_x:.03:max_x)
	ys = collect(min_y:.03:max_y)
	g_as_f = (a, b) -> g(x=>a, y=>b)
	zs = g_as_f.(xs', ys)
	PyPlot.contour(xs, ys, zs, levels=[0])
end


# ╔═╡ 0c6d7cdc-f8b4-11ea-14a6-9911030aba66
function plot_obstacle_setup(edge_size, a, b, eq_obstacles)
	PyPlot.clf()
	PyPlot.figure()
	PyPlot.xlim(-edge_size*1.1, edge_size*1.1)
	PyPlot.ylim(-edge_size*1.1, edge_size*1.1)
	for g=eq_obstacles
		plot_levelset(g, x...)
	end
	# PyPlot.plot([-edge_size, -edge_size, edge_size, edge_size, -edge_size],
	# 			[-edge_size, edge_size, edge_size, -edge_size, -edge_size],
	# 			ls="--", c="r")
	PyPlot.scatter([a[1],b[1]], [a[2], b[2]],  s=100, c="r")

end

# ╔═╡ 3f305e3e-f8b4-11ea-2ec4-8987fa78e852
begin
	offset = [0., 0]
	n = 2 # dimension of the space
	max_deg = 4
	num_pieces = 4
	reg = 0
	edge_size = 1.5
	a = [0.4, -1] .+ offset
	b = [0, 1].+ offset
	centers = [
		[-.4, .2],
		[.5, 0]
		]
	r = .5
end

# ╔═╡ 13e6d514-f8b4-11ea-3dc0-f7fe6a5c6f1c
begin
	@info "Model"

	model = Model(Ipopt.Optimizer)
	@variable(model, u[k=1:num_pieces, j=1:n], start=rand(num_pieces, n)[k,j])
	@variable(model, v[k=1:num_pieces, j=1:n], start=rand(num_pieces, n)[k,j])
	
	# Obstacle setup
	
	
end

# ╔═╡ 7a6efb50-f8b3-11ea-2616-c58441c1f6e2
begin
	
	
	
	
	
	
	#x = [ u[i,:] + t .* v[i,:] for i=1:num_pieces for t=0:.1:1]
	#circle_constraints = @expression model [i=1:size(x,1)] sum( (x[i] .-center).^2) - r^2
	lengths = @NLexpression model [i=1:size(v,1)] (sum(v[i, j]^2 for j=1:n))
	
	for center=centers
		for i=1:num_pieces
			for t=0:.1:1
				@NLconstraint model  sum( (u[i,j] + t*v[i,j]-center[j])^2 for j=1:n) >= r^2
			end
		end
	end
	
	for i=1:num_pieces-1
		@constraint model  u[i,:] .+ v[i,:] .== u[i+1,:]
	end
	
	@constraint model u[1,:] .== a
	@constraint model u[end,:] .+ v[end,:] .== b;
	
	
	for j=(1, -1)
	@constraint model [i=1:num_pieces] j .* u[i,:] .<= edge_size
	@constraint model [i=1:num_pieces] j .* u[i,:] .+ v[i, :] .<= edge_size
	end
	
	## Objective
	@info "Set Objective"
	@NLobjective model Min sum(l for l=lengths)
	@info "Optimize"
	optimize!(model)
	termination_status(model), objective_value(model)
	st = termination_status(model)
	st_opt = objective_value(model)
	
	
	@polyvar t
	opt_trajectory = [value.(Eμ.(u[i,:] .+ t .* v[i, :]) ) for i=1:num_pieces]
	@show(st, st_opt)
	@show(norm(a .-b))
	
	
end

# ╔═╡ 19570f66-f8b4-11ea-3715-439ae958d8b6
begin
	
	## Plot solution
	
	@polyvar x[1:n]
	plot_obstacle_setup(edge_size, a, b, [dist_squared(x, center) - r^2 for center=centers])
	for piece=opt_trajectory
		x1t = piece[1]
		x2t = piece[2]
		@show [x1t(0), x1t(1)], [x2t(0), x2t(1)]
		PyPlot.plot([x1t(0), x1t(1)], [x2t(0), x2t(1)])
	end
	PyPlot.title("NLP approach.")
	PyPlot.display_figs()
	
end

# ╔═╡ Cell order:
# ╠═89490f12-f8b3-11ea-267f-df907a3b7725
# ╟─06250748-f8b4-11ea-2dbc-492b649a4fa2
# ╟─08c04ee0-f8b4-11ea-0ae4-ebe7288fe0bf
# ╠═0c6d7cdc-f8b4-11ea-14a6-9911030aba66
# ╠═3f305e3e-f8b4-11ea-2ec4-8987fa78e852
# ╠═13e6d514-f8b4-11ea-3dc0-f7fe6a5c6f1c
# ╠═7a6efb50-f8b3-11ea-2616-c58441c1f6e2
# ╠═19570f66-f8b4-11ea-3715-439ae958d8b6
