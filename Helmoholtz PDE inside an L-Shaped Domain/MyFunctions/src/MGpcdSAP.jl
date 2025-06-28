function MGpcdSAP(
    K_global::SparseMatrixCSC{Float64, Int64},
    r_input::AbstractVector{Float64},  
    idx::Vector{Vector{Int}},
    mg_hier::Vector{MGHier},
    nPre::Int,
    nPost::Int,
    omega::Float64
    )

    nop = size(K_global,1);
    z = zeros(Float64, nop);

    for i = 1:3
        # Get subdomain rhs
        b_sub = r_input[idx[i]];

        z_sub = MG_Vcycle(b_sub, mg_hier[i], nPre, nPost, omega);
        z[idx[i]] .+= z_sub;
    end

    return z
end