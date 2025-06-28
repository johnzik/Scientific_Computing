function pcdSAP_vector_interface_wrapper(
    K_global_for_pcdSAP::SparseMatrixCSC{Float64, Int64}, 
    r_input_vec::AbstractVector{Float64},
    K_sub_all::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}}, 
    idx_all_subdomains::Vector{<:Vector{<:Integer}}
    )

    # Convert input 1D vector to an Nx1 matrix for the pcdSAP (originaly for pcdCG)
    r_input_matrix = reshape(r_input_vec, :, 1);
    
    # Call pcdSAP function
    z_output_matrix = pcdSAP(K_global_for_pcdSAP, r_input_matrix, K_sub_all, idx_all_subdomains);
    
    # Convert output Nx1 matrix back to a 1D vector
    return vec(z_output_matrix)
end