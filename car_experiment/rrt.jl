using LightGraphs

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

function rrt_find_path_tv(dim, start, goal, check_collision, check_collision_segment; num_iterations=10, step_size=.2, radius_goal=.1, obs_seed=0)

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
