function rrt_find_path(start, goal, check_collision, check_collision_segment; num_iterations=10, step_size=.2, radius_goal=.1)
	
	node_positions = [ start ]
	graph_edges = []
	Random.seed!(obs_seed)
	
	for counter=1:num_iterations
		
		# pick a random point
		new_x = 2 .* Random.rand(Float32, 2) .- 1
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
