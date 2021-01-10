using Ipopt

function shortest_path_nlp(
		n, num_cars, car_rad, obs_pos, obs_vel, obs_rad,
		a, b,
		num_pieces,
		edge_size,
		)
	@info "Model"


    model = Model(with_optimizer(Ipopt.Optimizer))
	@variable(model, u[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])
	@variable(model, v[k=1:num_pieces, j=1:n],
			  start=rand(num_pieces, n)[k,j])

	@info "Constraints"
		lengths = @NLexpression model [i=1:size(v,1)] sqrt(sum(v[i, j]^2 for j=1:n))

        # car avoidance
    for i=1:num_cars
        for j=(i+1):num_cars
            for s=1:num_pieces
                for t=0:.1:1
                    @NLconstraint model  (u[s, 2i-1] + t*v[s, 2i-1] - (u[s, 2j-1] + t*v[s, 2j-1]))^2 + (u[s, 2i] + t*v[s, 2i] - (u[s, 2j] + t*v[s, 2j]))^2  >= car_rad^2
                end
            end
        end
    end

    for i=1:num_cars
        for (op, ov)=zip(obs_pos, obs_vel)
            for s=1:num_pieces
                for t=0:.1:1
                    # time change variables
                    tt = (s-1 + t) / num_pieces
                    @NLconstraint model  (u[s, 2i-1] + t*v[s, 2i-1] - (op[1] + tt*ov[1]))^2 + (u[s, 2i] + t*v[s, 2i] - (op[2] + tt * ov[2]))^2  >= (car_rad + obs_rad)^2
                end
            end
        end
    end

    #
    # (t, x, y) -> (x-t)^2 + (y+t)^2 - obs_rad^2,
    # (t, x, y) -> (x+t)^2 + (y-t)^2 - obs_rad^2,
			# for i=1:num_pieces
			# 	for t=0:.1:1
			# 		obs_pos_t = obs_pos[j, :] .+ t .*  obs_vel[j, :]
			# 		@NLconstraint model  sum( (u[i,k] + t*v[i,k]-obs_pos_t[k] )^2 for k=1:n) >= radius_obs^2

			# 		@NLconstraint model  sum( (u[i,k] + t*v[i,k]-obs_pos_t[k] )^2 for k=1:n) >= radius_obs^2
			# 	end
			# end

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
