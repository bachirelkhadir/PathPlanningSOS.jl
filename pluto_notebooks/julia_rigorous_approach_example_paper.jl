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
	using DynamicPolynomials
	using PyPlot
	using LinearAlgebra
	using SumOfSquares
	using CSDP
	md"""A bunch of imports"""
end

# ╔═╡ c3d41f9e-f7de-11ea-2a56-b76588d6ef66
md"# Rigorous approach for predefined envs"

# ╔═╡ 5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
md"t $(@bind plot_at_time Slider(0:.01:1, show_value=true, default=0.))"

# ╔═╡ 98b2087e-02e4-11eb-31a3-1b9d90631910
md"# Problem data"

# ╔═╡ a511d81a-02e4-11eb-3419-4dfcb7b347a6
multiply_multipliers = false

# ╔═╡ 19278f2a-02e4-11eb-22fd-95719ebe4178
world_x = 1.01

# ╔═╡ 1b57adca-02e4-11eb-0821-d1159ff990dd
world_y = .995

# ╔═╡ 22ff0c08-02e4-11eb-08c6-ed077ca60192
start_x = -.99

# ╔═╡ 2b7f2fac-02e4-11eb-37cc-ebe8e56d810a
start_y = -.99

# ╔═╡ 2fa15284-02e4-11eb-0fe3-8179a9a62349
end_x = 0.012

# ╔═╡ 350361b8-02e4-11eb-0987-7965d90c3a30
end_y = .995

# ╔═╡ 3d5f7dbc-02e4-11eb-13af-8b7568f1f5af
num_pieces = 2

# ╔═╡ 448d1108-02e4-11eb-1390-c13a8f52358e
deg_relaxation_z = 10

# ╔═╡ 48d192ac-02e4-11eb-1f67-7b4f48aabc8a
log_reg = -8

# ╔═╡ 5339dc74-02e4-11eb-3e35-55a5370cedc8
obs_pos =  hcat([0., 0.25])'

# ╔═╡ 5e220858-02e4-11eb-1f4d-5b4f8e5b6ce5
obs_vel = hcat([0., .1])'

# ╔═╡ 66bd4bee-02e4-11eb-2852-ed82b19d422e
obs_radius = 5. /10.

# ╔═╡ 90161a4a-f7de-11ea-186d-972892ea3c26
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# ╔═╡ 708e3df4-02e4-11eb-3754-8750484484a5
	obstacles = [
    	(t, x) -> sum( (x .- obs_pos[i, :] .+ t .* obs_vel[i, :]).^2 ) - obs_radius^2 for i=1:size(obs_pos, 1)
	]

# ╔═╡ ab6bed40-fe8e-11ea-1787-8bfcbdf208c7
function find_path_using_rigorous_approach_cheap(n::Int, contraint_fcts, edge_size::Float64, 
    a::Array{Float64, 1}, b::Array{Float64, 1},
    max_deg_uv::Int, max_deg_z::Int, num_pieces::Int, solver,
	multiply_multipliers::Bool
    ;scale_init=1, reg=0, )
	# 	"""
	# 	Rigorous approach
	# 	No need for polyvar u since u = cumsum(v)
	# 	"""
    @polyvar v[1:num_pieces,1:n]
    @polyvar t
    @polyvar z[1:num_pieces]
    @polyvar x[1:n]
	u = vcat( (a .+ v[1,:] .* 0)',
			[a .+ sum([v[j, :] for j=1:i]) for i=1:num_pieces]'...)
	@info u
    @info "Definiting the model and the decision vars..."
    model = SOSModel(solver)
	@variable model γ
	
    contraint_polys = [
        f(t, x) for f in contraint_fcts
    ]

	
    @info "Declaration measures"
    μ, Eμ = PathPlanningSOS.measure_in_vars(model, [v..., z...], t, max_deg_uv, "μ")
    μz, Eμz = PathPlanningSOS.measure_in_vars(model, z, t, max_deg_z, "μz")
    decision_vars = [μ.a..., μz.a...]


    @info "Constraints definition"
    # total mass is one
    @constraint model Eμ(0*t+1).a[1] == 1
    @constraint model Eμz(0*t+1).a[1] == 1

    # obstacle constraints
    # localization of (u,v)
    loc_polynomials = [
	0*v[1] + 1,
	z...,
	#(sqrt(2)*edge_size .- z)...,
	[subs(g, [xj => xtj for (xj, xtj)=zip(x, u[i,:] .+ t .* v[i, :])]...)
			for i=1:num_pieces for g=contraint_polys]...
    ]

	if multiply_multipliers
    	loc_polynomials = [
     		f*g for f=loc_polynomials for g=loc_polynomials
	    ]
	end

    for g=loc_polynomials
	Mi_g = PathPlanningSOS.loc_matrix([v..., z...], t, g, max_deg_uv)
	M = Eμ.(Mi_g)
	var_M = hcat(map(Mij -> Mij.x, M)...)
        PathPlanningSOS.make_sos_on_0_1(model, t, M)
    end


    # x(0) ~ a, x(1) ~ b
	mons = monomials([v..., z...], 0:max_deg_uv-1)
    for m=mons
	# @constraint model Eμ.(m .* (u[1,:] .- a) ) .== 0
	@constraint model Eμ.(m .* (u[end,:].- b)) .== 0
    end


    # z^2 = ||v||^2
    mons = monomials([v..., z...], 0:max_deg_uv-2)
    for i=1:num_pieces
	@constraint model Eμ.(mons .* (z[i]^2-sum(v[i,:].^2))) .== 0
    end


    ## Constraint for μz

    for g=[0*z[1]+1, z...,]# (2*sqrt(2)*edge_size .- z)...]
	Mi_g = PathPlanningSOS.loc_matrix(z, t, g, max_deg_z)
	M = Eμz.(Mi_g)
	M = map(Mij -> Mij.a[1], M)
	if size(M,1) == 1
		@constraint model M[1] >= 0
	else
		@constraint model M in JuMP.PSDCone()
	end
	#PathPlanningSOS.make_psd(model, t, M)
    end


    # marginals of z agree
    mons_z = monomials(z, 0:max_deg_uv)
    @constraint model Eμ.(mons_z) .== Eμz.(mons_z)


    ## Objective
    @info "Set Objective"
    objective = sum(Eμ(zi) for zi=z)
    @constraint model [γ; decision_vars] in SecondOrderCone()
    objective = objective.a[1] + reg * γ
    @objective model Min objective

    #@objective model Min sum(Eμ(wi^4) for wi=[u..., v..., z...]).a[1]
    ##
    @info "Optimize"
    optimize!(model)
    @info termination_status(model), objective_value(model)
	
    opt_trajectory_pieces = [value.(Eμ.(u[i, :] + t .* v[i, :])) for i=1:num_pieces]
	opt_trajectory_length = sum([LinearAlgebra.norm(value.(Eμ.(v[i, :]))) 
		 						for i=1:num_pieces])
	
	@info "∑ E(z) = ", sum(value.(Eμ.(z)))
	@info "∑ ||E(v)|| = ",sum([LinearAlgebra.norm(value.(Eμ.(v[i, :]))) 
		 						for i=1:num_pieces])	
	@info "∑ √E(||v||^2) = ",sum([sqrt(value(Eμ.(sum(v[i, :].^2))).a[1]) 
		 						for i=1:num_pieces])	
	@info "∑ √E(z²) = ", sum(sqrt(value(Eμ(zi^2)).a[1]) for zi=z) 						
    return (objective_value(model),
			opt_trajectory_length,
			PathPlanningSOS.pieces_to_trajectory(opt_trajectory_pieces))
end

# ╔═╡ bd7b97e4-f7de-11ea-096f-27a885a176c7
begin
	result = []
	# compute optimal piece-wise linear trajectory
	for deg_relaxation=(2,3,4)
		n = 2
		reg = 10.0^log_reg
		opt_value, opt_trajectory_length, opt_trajectory = find_path_using_rigorous_approach_cheap(
			n, # n
			obstacles, # contraint_fcts
			world_x, # edge_size
			[start_x, start_y], # a
			[end_x, end_y], # b
			deg_relaxation, # max_deg_uv
			deg_relaxation_z, # max_deg_z
			num_pieces, 
			solver, 
			multiply_multipliers,
			reg=reg)
		@show deg_relaxation
		@show opt_value
		push!(result, deg_relaxation => opt_value)
	end
	result
end

# ╔═╡ 2e50f5d4-02e5-11eb-2655-c11875721100
begin
	deg_relaxation = 4
	n = 2
			reg = 10.0^log_reg
			opt_value, opt_trajectory_length, opt_trajectory = find_path_using_rigorous_approach_cheap(
				n, # n
				obstacles, # contraint_fcts
				world_x, # edge_size
				[start_x, start_y], # a
				[end_x, end_y], # b
				deg_relaxation, # max_deg_uv
				deg_relaxation_z, # max_deg_z
				num_pieces, 
				solver, 
				multiply_multipliers,
				reg=reg)
end

# ╔═╡ df5b5110-f7de-11ea-3a48-f15db9b1d873
begin
	q = PyPlot.figure()
	PathPlanningSOS.plot_at_time(plot_at_time, world_x, [start_x, start_y], [end_x, end_y],
						obstacles, opt_trajectory)
	PyPlot.title("t = $plot_at_time, relaxation deg $deg_relaxation, " *
				"($deg_relaxation_z en z), " *
				"mult-multipliers? $multiply_multipliers\n" *
				 "obj_value = $(round(opt_value, digits=2)), " *
				 "length = $(round(opt_trajectory_length, digits=2)), " *
	"reg = $reg")
	PyPlot.axes().set_aspect("equal")
	q
end

# ╔═╡ d6b5beea-02b2-11eb-23bb-d992de09a8b0
opt_value

# ╔═╡ df86eee0-02b2-11eb-1deb-6f9de8d2c31b
opt_trajectory_length

# ╔═╡ 6583ff0c-02ee-11eb-3fc9-ad85ee373074
opt_trajectory(.5)

# ╔═╡ caa28490-02e8-11eb-288a-4bf8121121a6
for i=(0,3,6,9)
	@show i, opt_trajectory(i/9)
end

# ╔═╡ Cell order:
# ╟─c3d41f9e-f7de-11ea-2a56-b76588d6ef66
# ╠═79575ad0-f7de-11ea-2e97-97d0955e0194
# ╠═5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
# ╟─df5b5110-f7de-11ea-3a48-f15db9b1d873
# ╠═d6b5beea-02b2-11eb-23bb-d992de09a8b0
# ╠═df86eee0-02b2-11eb-1deb-6f9de8d2c31b
# ╟─98b2087e-02e4-11eb-31a3-1b9d90631910
# ╟─a511d81a-02e4-11eb-3419-4dfcb7b347a6
# ╠═19278f2a-02e4-11eb-22fd-95719ebe4178
# ╟─1b57adca-02e4-11eb-0821-d1159ff990dd
# ╠═22ff0c08-02e4-11eb-08c6-ed077ca60192
# ╠═2b7f2fac-02e4-11eb-37cc-ebe8e56d810a
# ╠═2fa15284-02e4-11eb-0fe3-8179a9a62349
# ╠═350361b8-02e4-11eb-0987-7965d90c3a30
# ╠═3d5f7dbc-02e4-11eb-13af-8b7568f1f5af
# ╠═448d1108-02e4-11eb-1390-c13a8f52358e
# ╠═48d192ac-02e4-11eb-1f67-7b4f48aabc8a
# ╠═5339dc74-02e4-11eb-3e35-55a5370cedc8
# ╠═5e220858-02e4-11eb-1f4d-5b4f8e5b6ce5
# ╠═66bd4bee-02e4-11eb-2852-ed82b19d422e
# ╠═90161a4a-f7de-11ea-186d-972892ea3c26
# ╟─708e3df4-02e4-11eb-3754-8750484484a5
# ╟─ab6bed40-fe8e-11ea-1787-8bfcbdf208c7
# ╠═bd7b97e4-f7de-11ea-096f-27a885a176c7
# ╠═2e50f5d4-02e5-11eb-2655-c11875721100
# ╠═6583ff0c-02ee-11eb-3fc9-ad85ee373074
# ╠═caa28490-02e8-11eb-288a-4bf8121121a6
