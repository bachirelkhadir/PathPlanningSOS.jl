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


number pieces:
$(@bind num_pieces Slider(1:10, show_value=true, default=2))



relaxation degree:
$(@bind deg_relaxation Slider(1:8, show_value=true, default=2))

relaxation degree Z:
$(@bind deg_relaxation_z Slider(1:20, show_value=true, default=2))



Regularization:
$(@bind reg Slider(0.00:.01:1., show_value=true, default=.0))

Number obstacles:
$(@bind number_obs Slider(1:10, show_value=true, default=2))


Radius of obstacles:
$(@bind radius_obs Slider(0.00:.01:1, show_value=true, default=.3))

Speed obstacles
$(@bind speed_obs Slider(0.00:.1:1, show_value=true, default=0.))


obstacle seed:
$(@bind obs_seed Slider(1:100, show_value=true, default=0))


Multiply multipliers by themselves?
$(@bind multiply_multipliers CheckBox(default=false))


"""

# ╔═╡ 90161a4a-f7de-11ea-186d-972892ea3c26
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => false)

# ╔═╡ 96523c9a-f7de-11ea-1dd4-67eaad6f968d
begin
	Random.seed!(obs_seed)
	obs_pos = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obs_vel = 2 .* Random.rand(Float64, (number_obs, 2)) .- 1
	obstacles = [
    	(t, x) -> sum( (x .- obs_pos[i, :] .+ t .* speed_obs .* obs_vel[i, :]).^2 ) - radius_obs^2 for i=1:number_obs
	]
	
	md"Obstacle definitions here"
end

# ╔═╡ 2fe5cba8-fe86-11ea-083d-5d29545d976b
function find_path_using_rigorous_approach(n::Int, contraint_fcts, edge_size::Float64, 
    a::Array{Float64, 1}, b::Array{Float64, 1},
    max_deg_uv::Int, max_deg_z::Int, num_pieces::Int, solver,
    ;scale_init=1, reg=0, )
    @polyvar u[1:num_pieces,1:n]
    @polyvar v[1:num_pieces,1:n]
    @polyvar t
    @polyvar z[1:num_pieces]
    @polyvar x[1:n]

    @info "Definiting the model and the decision vars..."
    model = SOSModel(solver)
	@variable model γ
	
    contraint_polys = [
        f(t, x) for f in contraint_fcts
    ]

	
    @info "Declaration measures"
    μ, Eμ = PathPlanningSOS.measure_in_vars(model, [u..., v..., z...], t, max_deg_uv, "μ")
    μz, Eμz = PathPlanningSOS.measure_in_vars(model, z, t, max_deg_z, "μz")
    decision_vars = [μ.a..., μz.a...]


    @info "Constraints definition"
    # total mass is one
    @constraint model Eμ(0*t+1).a[1] == 1
    @constraint model Eμz(0*t+1).a[1] == 1

    # obstacle constraints
    # localization of (u,v)
    loc_polynomials = [
	0*u[1] + 1,
	z...,
	#(sqrt(2)*edge_size .- z)...,
	[subs(g, [xj => xtj for (xj, xtj)=zip(x, u[i,:] .+ t .* v[i, :])]...)
			for i=1:num_pieces for g=contraint_polys]...
    ]

    # loc_polynomials = [
    # 	f*g for f=loc_polynomials for g=loc_polynomials
    # ]

    for g=loc_polynomials
	Mi_g = PathPlanningSOS.loc_matrix([u..., v..., z...], t, g, max_deg_uv)
	M = Eμ.(Mi_g)
	var_M = hcat(map(Mij -> Mij.x, M)...)
        PathPlanningSOS.make_sos_on_0_1(model, t, M)
    end

    # continuity constraints
    # x(t=1+) = x(t=0-)
    mons = monomials([u..., v..., z...], 0:max_deg_uv-1)
    for i=1:(num_pieces-1)
	end_piece_i = u[i,:] .+ v[i, :]
	start_piece_iplus = u[i+1,:]
	for m=mons
 	    @constraint(model, Eμ.(m .* (end_piece_i .- start_piece_iplus)) .== 0)
	end
    end
    # x(0) ~ a, x(1) ~ b
    for m=mons
	@constraint model Eμ.(m .* (u[1,:] .- a) ) .== 0
	@constraint model Eμ.(m .* (u[end,:]+v[end,:] .- b)) .== 0
    end


    # z^2 = u^2+v^2
    mons = monomials([u..., v..., z...], 0:max_deg_uv-2)
    for i=1:num_pieces
	@constraint model Eμ.(mons .* (z[i]^2-sum(v[i,:].^2))) .== 0
    end


    ## Constraint for μz


    # psd
    # Mz = loc_matrix(z, 0*z[1]+1, max_deg_z)
    # Mz = Eμz.(Mz)
    # Mz = map(Mij -> Mij.a[1], Mz)
    # @constraint model Mz in PSDCone()
    #

    for g=[0*z[1]+1, z..., (sqrt(2)*edge_size .- z)...]
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
    PathPlanningSOS.pieces_to_trajectory(opt_trajectory_pieces)
end


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
	@info "∑ √E(z²) = ", sum(sqrt(value(Eμ(zi^2)).a[1]) for zi=z) 						
    return (objective_value(model),
			opt_trajectory_length,
			PathPlanningSOS.pieces_to_trajectory(opt_trajectory_pieces))
end

# ╔═╡ bd7b97e4-f7de-11ea-096f-27a885a176c7
begin
	# compute optimal piece-wise linear trajectory
	n = 2
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
	md"Computing Optimal path..."
end

# ╔═╡ 5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
md"t $(@bind plot_at_time Slider(0:.01:1, show_value=true, default=0.))"

# ╔═╡ df5b5110-f7de-11ea-3a48-f15db9b1d873
begin
	q = PyPlot.figure()
	PathPlanningSOS.plot_at_time(plot_at_time, world_x, [start_x, start_y], [end_x, end_y],
						obstacles, opt_trajectory)
	PyPlot.title("t = $plot_at_time, relaxation deg $deg_relaxation, " *
				"($deg_relaxation_z en z), " *
				"mult-multipliers? $multiply_multipliers\n" *
				 "obj_value = $(round(opt_value, digits=2)), " *
				 "length = $(round(opt_trajectory_length, digits=2))")
	PyPlot.axes().set_aspect("equal")
	q
end

# ╔═╡ Cell order:
# ╟─c3d41f9e-f7de-11ea-2a56-b76588d6ef66
# ╟─79575ad0-f7de-11ea-2e97-97d0955e0194
# ╟─df5b5110-f7de-11ea-3a48-f15db9b1d873
# ╟─8560dd74-f7de-11ea-3c51-7f36d640d43b
# ╟─90161a4a-f7de-11ea-186d-972892ea3c26
# ╟─96523c9a-f7de-11ea-1dd4-67eaad6f968d
# ╟─2fe5cba8-fe86-11ea-083d-5d29545d976b
# ╟─ab6bed40-fe8e-11ea-1787-8bfcbdf208c7
# ╠═bd7b97e4-f7de-11ea-096f-27a885a176c7
# ╟─5e97a0fa-f7df-11ea-1750-d7ce2d803e9d
