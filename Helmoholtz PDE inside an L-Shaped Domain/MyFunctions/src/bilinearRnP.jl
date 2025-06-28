"""
    bracketing1D(val::Float64, sorted_lines::Vector{Float64})

This is a helper function made to create the 1D bracketing, i.e. to find the lines that bracket the fine node and
output the Normalized Interpolation Coord xi

# Arguments
- **val::Float64 :** The fine node coordinate (x or y axis)
- **sorted_lines::Vector{Float64} :** The sorted coarse x or y lines

# Returns
- **idx1, idx2 :** The 2 coarse lines that bracket the fine node
- **xi :** The Normalized Interpolation Coordinate
"""
function bracketing1D(val::Float64, sorted_lines::Vector{Float64})
    len_lines = length(sorted_lines)
    if len_lines == 0
        return 0, 0, 0.0 # No lines to define a segment
    end
    if len_lines == 1 # Only one coarse line available
        # Return 1, 1 bracketing coarse lines &
        # Check if fine node (val) ≈ coordinate of the single coarse line -> return 0.0 (normalized pos)
        # Else check if val < coarse line's position -> if yes return 0.0, if not return 1.0
        return 1, 1, (isapprox(val, sorted_lines[1]) ? 0.0 : (val < sorted_lines[1] ? 0.0 : 1.0))
    end

    # Find k such that sorted_lines[k-1] <= val <= sorted_lines[k] (Bracketing)
    k_insert = searchsortedfirst(sorted_lines, val);

    idx1, idx2 = 0, 0; # Initialize 1st and 2nd bracketing coarse line indexes, respectively

    if k_insert == 1 # val <= sorted_lines[1] (at or before the very first coarse line)
        idx1 = 1;
        idx2 = 2; # Use the first segment [lines[1], lines[2]] for interpolation

    elseif k_insert > len_lines # val > sorted_lines[end] (after the very last coarse line)

        idx1 = len_lines - 1;
        idx2 = len_lines; # Use the last segment [lines[end-1], lines[end]] for interpolation

    else # sorted_lines[k_insert-1] < val <= sorted_lines[k_insert]

        if isapprox(val, sorted_lines[k_insert]) # val is on line[k_insert]
            idx1 = k_insert;
            idx2 = k_insert; # Coincident with this line
        else # val is between line[k_insert-1] and line[k_insert]
            idx1 = k_insert - 1;
            idx2 = k_insert;
        end
    end
    
    # If coincident, xi depends on which end of segment we consider it
    # If fine node is on coarse_lines[idx1], weights involve xi=0
    # If fine node is on coarse_lines[idx2], weights involve xi=1
    if idx1 == idx2 # Coincident with sorted_lines[idx1]
        # To form a segment for interpolation, if it's not the last point, use [idx1, idx1+1] with xi=0
        # If it is the last point, use [idx1-1, idx1] with xi=1
        if isapprox(val, sorted_lines[idx1])
            if idx1 < len_lines
                return idx1, idx1 + 1, 0.0 # Treat as start of segment [idx1, idx1+1]
            else
                return idx1 - 1, idx1, 1.0 # Treat as end of segment [idx1-1, idx1]
            end
        else # Should not happen if idx1==idx2 was due to coincidence (Error)
             return 0, 0, 0.0
        end
    end
    
    line_val1 = sorted_lines[idx1];
    line_val2 = sorted_lines[idx2];
    denominator = line_val2 - line_val1;

    if abs(denominator) < 1e-12 # Avoid division by zero
        # If val is also on this line, treat xi as 0 or 1 depending on which index it matched more closely
        xi_norm = isapprox(val, line_val1) ? 0.0 : (isapprox(val, line_val2) ? 1.0 : 0.5) # Fallback if lines too close
    else
        xi_norm = (val - line_val1) / denominator
    end
    
    return idx1, idx2, clamp(xi_norm, 0.0, 1.0) # Clamp to handle extrapolation robustly
end

"""
    build_P_bilinear_and_R_transpose(fine_coords::Matrix{Float64})

This function builds and returns the Prolongation and the Restriction operators and the list of fine 
node indices that are also coarse nodes.

# Arguments
- **fine_coords::Matrix{Float64} :** All the fine (initial) coordinates of the grid

# Returns
- **R_op :** The Restriction operator
- **P_op :** The Prolongation operator
- **fine_indices_selected_as_coarse :** The list of fine node indices that are also coarse nodes
"""
function build_P_bilinear_and_R_transpose(fine_coords::Matrix{Float64})
    n_fine = size(fine_coords, 1);

    minNodes = 4 # Min fine nodes to attempt meaningful 2x2 coarse cell
    if n_fine < minNodes
        @warn "P-bilinear: Grid too small ($n_fine nodes). Returning identity operators."
        identity_op = sparse(1.0I, n_fine, n_fine);
        return identity_op, identity_op', collect(1:n_fine)
    end

    unique_x_fine = sort(unique(fine_coords[:, 1]));
    unique_y_fine = sort(unique(fine_coords[:, 2]));

    if length(unique_x_fine) < 2 || length(unique_y_fine) < 2
        @warn "P-bilinear: Not enough unique fine grid lines for coarsening. Returning identity."
        identity_op = sparse(1.0I, n_fine, n_fine);
        return identity_op, identity_op', collect(1:n_fine)
    end

    # Define Coarse Grid structure from fine_coords
    coarse_x_lines = unique_x_fine[1:2:end];
    coarse_y_lines = unique_y_fine[1:2:end];

    if length(coarse_x_lines) < 1 || length(coarse_y_lines) < 1 # Need at least one line
        @warn "P-bilinear: Coarsening produced no coarse grid lines. Returning identity."
        identity_op = sparse(1.0I, n_fine, n_fine);
        return identity_op, identity_op', collect(1:n_fine)
    end
    
    # Create a map from coarse (x,y) to its 1D local coarse index (1 to n_coarse)
    # And identify the original fine_coords indices that are these coarse nodes
    map_coarse_coord_to_local_idx = Dict{Tuple{Float64, Float64}, Int}();
    fine_indices_selected_as_coarse = Int[]; # Will store original fine indices
    

    temp_coarseCoords = Tuple{Float64,Float64}[]; # Coordinates of actual coarse nodes
    map_fine_coord_to_fine_idx = Dict((fine_coords[i,1], fine_coords[i,2]) => i for i=1:n_fine); # Fine coords to idx mapping

    for y_val in coarse_y_lines
        for x_val in coarse_x_lines
            # Find fine nodes that are these coarse_x_lines/coarse_y_lines intersections
             fine_idx_match = get(map_fine_coord_to_fine_idx, (x_val, y_val), 0);
             if fine_idx_match != 0
                push!(temp_coarseCoords, (x_val,y_val));
             end
        end
    end

    # Re-sort temp_coarseCoords to ensure canonical order for map_coarse_coord_to_local_idx
    sort!(temp_coarseCoords, by = p->(p[2],p[1]));
    
    n_coarse = length(temp_coarseCoords);
    if n_coarse == 0 || n_coarse == n_fine # No actual coarse nodes formed or no coarsening
         @warn "P-bilinear: No effective coarse grid defined or no coarsening. Returning identity.";
        identity_op = sparse(1.0I, n_fine, n_fine);
        return identity_op, identity_op', collect(1:n_fine)
    end

    for i=1:n_coarse
        map_coarse_coord_to_local_idx[temp_coarseCoords[i]] = i;
        # Populate fine_indices_selected_as_coarse based on this sorted order
        push!(fine_indices_selected_as_coarse, map_fine_coord_to_fine_idx[temp_coarseCoords[i]]);
    end

    P_I, P_J, P_V = Int[], Int[], Float64[]; # Triplets for sparse P

    for i_fine = 1:n_fine
        xf, yf = fine_coords[i_fine, 1], fine_coords[i_fine, 2];

        ix1_line_idx, ix2_line_idx, xi_norm = bracketing1D(xf, coarse_x_lines);
        iy1_line_idx, iy2_line_idx, yi_norm = bracketing1D(yf, coarse_y_lines);
        
        # Define the 4 parent coarse cell corner coordinates
        parent_CornerCoords = [
            (coarse_x_lines[ix1_line_idx], coarse_y_lines[iy1_line_idx]), # c00
            (coarse_x_lines[ix2_line_idx], coarse_y_lines[iy1_line_idx]), # c10
            (coarse_x_lines[ix1_line_idx], coarse_y_lines[iy2_line_idx]), # c01
            (coarse_x_lines[ix2_line_idx], coarse_y_lines[iy2_line_idx])  # c11
        ];

        weights = [
            (1 - xi_norm) * (1 - yi_norm),  # for c00
            xi_norm       * (1 - yi_norm),  # for c10
            (1 - xi_norm) * yi_norm,        # for c01
            xi_norm       * yi_norm         # for c11
        ];
        
        total_weight = 0.0;
        active_parents_data = Tuple{Int, Float64}[]; # (coarse_local_idx, weight)

        for k_corner = 1:4
            coarseNodes = parent_CornerCoords[k_corner];
            weight = weights[k_corner];
            
            coarse_node_local_idx = get(map_coarse_coord_to_local_idx, coarseNodes, 0);
            if coarse_node_local_idx != 0 && weight > 1e-9 # If parent coarse node exists in the list and weight is significant
                push!(active_parents_data, (coarse_node_local_idx, weight));
                total_weight += weight;
            end
        end
        
        if total_weight < 1e-9 # If no valid parents or all weights effectively zero
            min_d_sq = Inf; # Smallest squared distance so far
            best_c_local_idx = 0; # Stores the idx of the closest coarse node
            for i_c_search = 1:n_coarse
                c_coord = temp_coarseCoords[i_c_search];
                d_sq = (xf-c_coord[1])^2 + (yf-c_coord[2])^2; # Euclidean distance (squared distance)
                if d_sq < min_d_sq 
                    min_d_sq = d_sq; 
                    best_c_local_idx = i_c_search; 
                end
            end
            if best_c_local_idx > 0
                push!(P_I, i_fine); push!(P_J, best_c_local_idx); push!(P_V, 1.0);
            end
            continue
        end

        # Normalize weights for this fine node to sum to 1
        for (coarse_idx, weight) in active_parents_data
            push!(P_I, i_fine);
            push!(P_J, coarse_idx); # This is the local coarse index (1 to n_coarse)
            push!(P_V, weight / total_weight) ;
        end
    end

    if isempty(P_I)
        @warn "P-bilinear: Prolongation operator P is empty. Returning identity.";
        identity_op = sparse(1.0I, n_fine, n_fine);
        return identity_op, identity_op', collect(1:n_fine)
    end
    
    P_op = sparse(P_I, P_J, P_V, n_fine, n_coarse);
    R_op = P_op'; # R = P^T

    return R_op, P_op, fine_indices_selected_as_coarse
end