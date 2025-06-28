"""
    SubdomainOverlap_Matrices(GO::Vector{Matrix{Int64}}, K::SparseMatrixCSC{Float64, Int64} , F_global::Matrix{Float64}, H::Int)

This function takes the Global Ordering, in subdomains, and overlaps each into the other by H elements. Then 
it constructs and returns the stiffness and right hand side matrices of each subdomain.

# Arguments
- **GO::Vector{Matrix{Int64}} :** The Global Ordering matrix which contains the 3 seperate subdomains matrices meshed (use Lmesh to create it)
- **K::SparseMatrixCSC{Float64, Int64} :** K is the sparse stiffness matrix of our inhomogeneous Helmholtz equation
- **F_global::Matrix{Float64} :** F_global is the final right hand side of the equation. Now use a solver to solve the system K*u=F_global for u(x,y,t)
- **H::Int :** This variable sets the overlap number, in elements, between the subdomains

# Returns
- **K_sub :** This contains the stiffness matrices for each subdomain (outputs: K_sub[i], i=1,2,3)
- **F_global_sub :** This contains the right hand side matrices for each subdomain (outputs: F_global_sub[i], i=1,2,3)
- **idx :** The 3 subdomains node vectors with their overlaps
"""
function SubdomainOverlap_Matrices(GO::Vector{Matrix{Int64}}, 
                                   K::SparseMatrixCSC{Float64, Int64}, 
                                   F_global::Matrix{Float64}, 
                                   H::Int)

    # Catch Error: Overlap must not exceed the element length of a subdomains vertical or horizontal side
    if H > length(GO[1][:,1]) - 1
        return error("Overlap H must NOT exceed the element length of a subdomains vertical or horizontal side");
    end

    # Inner node vectors
    inGO1 = vec(GO[1]);
    inGO2 = vec(GO[2]);
    inGO3 = vec(GO[3]);

    # Overlapping node vectors
    G12 = vec(GO[2][:, 1:H+1]); # Overlap of subdomain 1 into 2 (H+1 because H is the overlap in elements)
    G13 = vec(GO[3][end-H:end, :]); # Overlap of subdomain 1 into 3
    G1 = vcat(G12, G13); # Merge the overlapping nodes of subdomain 1 into 2 & 3
    G21 = vec(GO[1][:, end-H:end]); # Overlap of subdomain 2 into 1
    G31 = vec(GO[1][1:H+1, :]); # Overlap of subdomain 3 into 1

    # Overlapping and inner node union vectors
    idx1 = unique([inGO1;G1]);
    idx2 = unique([inGO2;G21]);
    idx3 = unique([inGO3;G31]);

    # Subdomain stiffness matrices
    K_sub = [
        K[idx1, idx1],
        K[idx2, idx2],
        K[idx3, idx3]
    ]

    # Subdomain right hand side
    F_global_sub = [
        F_global[idx1],
        F_global[idx2],
        F_global[idx3]
    ]

    return K_sub, F_global_sub, [idx1, idx2, idx3]
end