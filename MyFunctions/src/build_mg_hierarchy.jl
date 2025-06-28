function build_mg_hierarchy(
    A_fine::SparseMatrixCSC{Float64, Int64}, 
    coords_arg::Matrix{Float64},
    nodes_arg::Vector{Int},
    min_size::Int=5, 
    max_levels::Int=5
    )

    levels = MGLevel[];
    current_A = A_fine;
    current_coords = coords_arg;
    current_nodes = nodes_arg;
    
    for _ in 1:max_levels

        if size(current_A, 1) <= min_size
            println("Current A size ($(size(current_A,1))) <= min_size ($min_size). Stopping coarsening.");
            break
        end
        
        # Build injection Restriction & Prolongation operators
        R, P, coarse_indices = build_P_bilinear_and_R_transpose(current_coords);
        
        nc = size(R, 1) # Number of coarse nodes for this new level
        
        if nc <= min_size
            println("Number of coarse nodes ($nc) <= min_size ($min_size) after injection. Stopping coarsening.");
            break
        end
        
        # Galerkins coarse operator
        A_coarse = R * current_A * P;
        
        # Store current level's A and R, P
        push!(levels, MGLevel(current_A, R, P, current_nodes, current_coords));
        
        # Prepare next level
        current_A = A_coarse;
        current_coords = current_coords[coarse_indices, :];
        current_nodes = current_nodes[coarse_indices];
    end
    
    # Add coarsest level
    push!(levels, MGLevel(current_A, 
                         nothing,  # R = nothing
                         nothing,  # P = nothing
                         current_nodes,
                         current_coords));
    
    return MGHier(levels)
end
