module PathPlanningSOS


export find_path_using_heuristic
export show_animation_in_notebook
export make_animation_as_video
export make_animation_as_plots


using Base64
using DynamicPolynomials
using LinearAlgebra
using MultivariateMoments
using Plots
using ProgressMeter
using PyCall
using PyPlot
using Random
using SumOfSquares
using Test
pyanim = PyPlot.matplotlib.animation




## Helper functions
dist_squared(a, b) = sum((a .- b).^2)


"""
Returns the localization matrix (g * mi * mj)_ij,
where mi and mj are monomials in `vars` up to degree `max_deg`.
"""

function loc_matrix(vars, t, g, max_deg)
    half_deg = (max_deg - maxdegree(subs(g, t=>1))) ÷ 2
    half_mons = monomials(vars, 0:half_deg)
    M = g .* (half_mons * half_mons')
    M
end


"""
Returns a newly created measure `ϕ` in variables `vars`
"""

function measure_in_vars(model, vars, t, max_deg, name)
    mons = monomials(vars, 0:max_deg)
    moments = @variable(model, [i=1:size(mons,1)], base_name=name)
    μ = MultivariateMoments.measure(moments, mons)
    E_μ = p -> MultivariateMoments.expectation(μ, p)
    # make E_μ handle polynomials that depend on t
    E_μ_tv = p -> begin
        mons_p = monomials(p)
        coeff_p = coefficients(p)
        r = 0
        for (m, c)=zip(mons_p, coeff_p)
            deg_t = degree(m, t)
            m_without_t = subs(m, t=>1)
            r += c * E_μ(m_without_t) * t^deg_t
        end
        r
    end

    return μ, E_μ_tv
end

"""
Make the matrix M(t) sos on [0, 1].
"""
function make_sos_on_0_1(model, t, M)
    # polynomial variables inside M
    var_M = hcat(map(Mij -> Mij.x, M)...)

    # if M does not depent on `t`,
    # there is no need to invoke the SOS machinery
    if t in var_M
        y = [similarvariable(eltype(M), gensym()) for i in 1:size(M, 1)]
        p = dot(y, M * y)
        # todo: adapt multipliers to degree of t in M
        domain_t = @set t*(1-t) >= 0
        @constraint(model, p >= 0, domain=domain_t )
    else
        # hack to convert a constant polynomial to a constant
        M = map(Mij -> Mij.a[1], M)
        if size(M,1) == 1
            @constraint model M[1] >= 0
        else
            @constraint model M in JuMP.PSDCone()
        end
    end
end


function find_path_using_heuristic(n, contraint_fcts, edge_size, a, b,
    max_deg_uv, num_pieces, solver,
    weight_lenght,
    num_iterations,
    ;scale_init=1, reg=0, seed=0)


    @polyvar u[1:n]
    @polyvar v[1:n]
    @polyvar t
    @polyvar x[1:n]
    xt = u .+ t .* v
    measure_vars = [u..., v...]
    contraint_polys = [
        f(t, x) for f in contraint_fcts
    ]


    @info "Definiting the model and the decision vars..."
    model = SOSModel(solver)

    @variable model γ # γ = ||decision_vars||
    @variable model α[1:num_pieces] # α = ||E(v)||

    μ_uvs, Eμ_uvs = zip([
        measure_in_vars(model, measure_vars, t, max_deg_uv, "uv_$i")
        for i=1:num_pieces
            ]...)

    decision_vars = cat([μ.a for μ=μ_uvs]..., dims=1)

    ## Constraints
    @info "Constraints definition..."


    # total mass is one
    for Eμ=Eμ_uvs
        @constraint model Eμ(0*t+1) == 1
    end


    # obstacle constraints
    # localization of (u,v)

    loc_polynomials = [
        0*u[1] + 1,
        contraint_polys...
    ]

    for (i, E_uv)=enumerate(Eμ_uvs)

        for g=loc_polynomials
            g_time_corrected = subs(g, t=>(t+i-1)/num_pieces,
                                [xj => xtj for (xj, xtj)=zip(x, xt)]...)


            Mi_g = loc_matrix(measure_vars, t, g_time_corrected, max_deg_uv)
            M = E_uv.(Mi_g)

            # make sos on [0, 1]
            make_sos_on_0_1(model, t, M)

        end
    end


    # continuity constraints
    # x(t=1+) = x(t=0-)

    for (E_uv, E_uv_plus)=zip(Eμ_uvs, Eμ_uvs[2:end])
        c = @constraint(model, E_uv.(subs(xt, t=>1)) .== E_uv_plus.(subs(xt, t=>0)))
    end



    # x(0) ~ a, x(1) ~ b
    mons = monomials([u..., v...], 0:max_deg_uv-1)
    for m=mons
        @constraint model (Eμ_uvs[1]).(m .* (subs(xt, t=>0) .- a)) .== 0
        @constraint model (Eμ_uvs[end]).(m .* (subs(xt, t=>1) .- b)) .== 0
    end

    # regularization

    @constraint model [γ; decision_vars] in SecondOrderCone()

    for i=1:num_pieces
        piece_i = [m.a[1] for m=Eμ_uvs[i].(v)]
        @constraint model [α[i]; piece_i...] in SecondOrderCone()
    end

    ## Random initialization


    Random.seed!(seed)
    uv_k = []

    push!(uv_k, scale_init .* (2 .* Random.rand(Float64, (num_pieces, n)) .- 1))

    @info "Starting iterative heuristic..."
    @showprogress for k=1:num_iterations
        objective = sum([ Eμ(ui^2) - 2*Eμ(ui)*old_ui
                for (Eμ,old_u)=zip(Eμ_uvs, uv_k[end])
                for (old_ui,ui)=zip(old_u, measure_vars)])

        objective = objective.a[1] + reg * γ + weight_lenght * sum(α)
        @objective model Min objective

        optimize!(model)
        termination_status(model), objective_value(model)
        st = termination_status(model)
        st_opt = objective_value(model)
        opt_trajectory = [value.(Euv.(xt)) for Euv in Eμ_uvs]

        push!(uv_k,  [ [value(Eμ(ui))
                        for ui=measure_vars] for Eμ=Eμ_uvs])

        end
        opt_trajectory_pieces = [value.(Euv.(xt)) for Euv in Eμ_uvs]
        pieces_to_trajectory(opt_trajectory_pieces)
end


"""
Convert linear pieces to a single trajectory.
"""

function pieces_to_trajectory(trajectory_pieces)
     t -> begin
        num_pieces = size(trajectory_pieces, 1)
        idx_piece = Int(fld(t * num_pieces, 1))
        idx_piece = min(idx_piece, num_pieces-1)
        piece = trajectory_pieces[idx_piece+1]
        # time shift
        map(xi -> xi(t*num_pieces-idx_piece), piece)
    end
end

# Animation helpers

function plot_levelset(g)
    """
    Plot the set {x | g(x) = 0}
    """
    min_x, max_x = PyPlot.xlim()
    min_y, max_y = PyPlot.ylim()
    xs = collect(min_x:.03:max_x)
    ys = collect(min_y:.03:max_y)
    zs = g.(xs', ys)
    PyPlot.contour(xs, ys, zs, levels=[0])
end


function plot_at_time(t, edge_size, a, b, eq_obstacles, opt_trajectory)
    PyPlot.xlim(-edge_size*1.1, edge_size*1.1)
    PyPlot.ylim(-edge_size*1.1, edge_size*1.1)

    ts = collect(0:0.01:1)
    xts = hcat(opt_trajectory.(ts)...)
    xt = opt_trajectory(t)
    PyPlot.plot(xts[1, :], xts[2, :])
    PyPlot.scatter([xt[1]], [xt[2]], s=100, c="g")


    for g=eq_obstacles
        plot_levelset((x, y) -> g(t, [x, y]))
    end

    PyPlot.scatter([a[1],b[1]], [a[2], b[2]],  s=100, c="r")
end


function show_animation_in_notebook(filename)
    base64_video = base64encode(open(filename))
    display("text/html", """<video controls src="data:video/x-m4v;base64,$base64_video">""")
end


function make_animation_as_video(filename, opt_trajectory,
        constraint_fcts, edge_size, a, b;
        num_frames=21, interval_between_frames=100)

    @info "Storing animation in $filename ..."

    fig = PyPlot.figure(figsize=(5,5))


    function make_frame(i)
        fig.clf()
        t = i/num_frames
        plot_at_time(t, edge_size, a, b, constraint_fcts, opt_trajectory)
        PyPlot.title("t = $t")
    end

    PyPlot.withfig(fig) do
        myanim = pyanim.FuncAnimation(fig, make_frame,
            frames=num_frames, interval=interval_between_frames)
        myanim.save(filename, bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    end
end

function make_animation_as_plots(opt_trajectory,
        constraint_fcts, edge_size, a, b;
        num_frames=21)

    @info "Starting plots ..."

    fig = PyPlot.figure(figsize=(5,5))
    @showprogress for i=0:num_frames
        #fig.cla()
        t = i/num_frames
        plot_at_time(t, edge_size, a, b, constraint_fcts, opt_trajectory)
        PyPlot.title("t = $t")
        fig.clear()
        fig.show()
    end
end

end # module
