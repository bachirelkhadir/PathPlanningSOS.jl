
# opt_trajectory is a function that takes t as input,
# and gives the location of the particle at time t.
#
function plot_sequence(opt_trajectory, edge_size, g_obstacles)
    fig = figure("Car movement",figsize=(20,2)) # Create a new blank figure
    N = 8
    car_colors = ["r", "g", "b", "m"]

    full_trajectory = hcat(opt_trajectory.(0:.05:1)...)

    for t=0:(N-1)
        subplot(101 + N*10 + t, aspect="equal")
        PyPlot.title("t = $(t/(N-1))")
        PyPlot.xlim(-edge_size*1.1, edge_size*1.1)
        PyPlot.ylim(-edge_size*1.1, edge_size*1.1)
        # PyPlot.axis("equal")
        xt = opt_trajectory.(1/(N-1) * t)

        # plt.gcf().gca().add_artist(plt.Circle((0, 0), radius_obs, fill=true))
        for g ∈ g_obstacles
            PathPlanningSOS.plot_levelset((x, y) -> g(1/(N-1) * t, x, y))
        end
        for i=1:num_cars
            plt.gcf().gca().add_artist(plt.Circle(xt[2*i-1:2*i], radius_car, color=car_colors[i], fill=true))
            if t == 0
                PyPlot.plot(full_trajectory[2*i-1, :], full_trajectory[2*i, :], alpha=.5, color=car_colors[i], ls="--")
            end
        end
    end
end


function plot_animation(opt_trajectory, edge_size, g_obstacles)

    car_colors = ["r", "g", "b", "m"]
    @info "Generating animation"
    @showprogress for (i,t)=enumerate(0:.05:1)
        PyPlot.figure()
        PyPlot.axis("equal")
        PyPlot.title("t = $t")

        PyPlot.xlim(-edge_size*1.1, edge_size*1.1)
        PyPlot.ylim(-edge_size*1.1, edge_size*1.1)
        xt = opt_trajectory.(t)


        for g=g_obstacles
        PathPlanningSOS.plot_levelset((x, y) -> g(t, x, y))
        end
        for i=1:num_cars
            plt.gcf().gca().add_artist(plt.Circle(xt[2*i-1:2*i], radius_car, color=car_colors[i], fill=true))
        end

        fn_index = lpad(i, 3, "0")
        PyPlot.savefig("imgs/path_planning_frame_$(fn_index).png")
    end
end
