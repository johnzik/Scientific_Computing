function MGHier_Setup(
    K_sub::Vector{SparseMatrixCSC{Float64, Int64}}, 
    idx::Vector{Vector{Int64}}, 
    all_coords::Matrix{Float64},
    num_mg_levels::Int,
    min_coarse_size_param::Int
    )
    println("\n--- Setting up Multigrid Hierarchies for All Subdomains ---")
    mg_hierarchies_list = Vector{MGHier}(undef, length(K_sub)) # To store MGHier for each subdomain

    for i = 1:3
        println("  Processing Subdomain $i for MG setup...");
        A_sub_i = K_sub[i];
        subdomain_global_indices = idx[i];

        if isempty(A_sub_i) || isempty(subdomain_global_indices)
            @warn "  Subdomain $i: A_sub or node list is empty. Cannot build hierarchy."
            # Create an MGHier with no levels, or a single level if A_sub_i exists
            mg_hierarchies_list[i] = MGHier(MGLevel[]) # Empty hierarchy
            continue
        end

        # Get coordinates for the current subdomain's nodes
        sub_coords_i = all_coords[subdomain_global_indices, :]

        mg_hierarchies_list[i] = build_mg_hierarchy(
            A_sub_i,
            sub_coords_i,
            copy(subdomain_global_indices), # Pass a copy of original global indices for this subdomain
            min_coarse_size_param,
            num_mg_levels
        );

        if isempty(mg_hierarchies_list[i].levels)
            @warn "  Subdomain $i: MG hierarchy is empty after build."
        else
            println("  Subdomain $i: Stored $(length(mg_hierarchies_list[i].levels)) MG levels.")
        end
    end
    println("--- Finished Multigrid Hierarchies Setup ---");
    return mg_hierarchies_list;
end