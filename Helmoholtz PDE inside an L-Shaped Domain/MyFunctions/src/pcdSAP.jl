"""
    pcdSAP(K::SparseMatrixCSC{Float64, Int64}, 
           r_input::Matrix{Float64}, 
           K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}}, 
           idx::Vector{Vector{Int}})

This function is a preconditioner, utilizing Schwartz Additive Procedure, that is used
together with the preconditioned Conjugate Gradient Method (pcdCG function) and GMRes (pcdGMRes function)

# Arguments
- **K::SparseMatrixCSC{Float64, Int64} :** K is the sparse stiffness matrix of our inhomogeneous Helmholtz equation
- **r_input::Matrix{Float64} :** This is the residual input vector of the pcdCG function 
- **K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}} :** This contains the stiffness matrices for each subdomain factorized (outputs: K_sub[i], i=1,2,3)
- **idx::Vector{Vector{Int}} :** The 3 subdomains node vectors with their overlaps

# Returns
 - **z :** This is the solution of the system z = M^-1 * r , where M^-1 is the preconditioner matrix
"""
function pcdSAP(K::SparseMatrixCSC{Float64, Int64}, 
                r_input::Matrix{Float64}, 
                K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}}, 
                idx::Vector{Vector{Int}} 
                )

    nop = size(K,1); # Number of nodes in the whole domain

    z = zeros(nop, 1); # Output vector initialization

    # Solve system for the union of subdomain 1 and its overlaps
    z1 = K_sub[1] \ r_input[idx[1], 1];

    # Solve system for the union of subdomain 2 and its overlaps
    z2 = K_sub[2] \ r_input[idx[2], 1];

    # Solve system for the union of subdomain 3 and its overlaps
    z3 = K_sub[3] \ r_input[idx[3], 1];

    z[idx[1], 1] .+= z1;
    z[idx[2], 1] .+= z2;
    z[idx[3], 1] .+= z3;

    return z;
end