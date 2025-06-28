"""
    pcdCG(K::SparseMatrixCSC{Float64, Int64}, 
               F::Matrix{Float64}, 
               tol::Float64, 
               K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}},
               idx::Vector{Vector{Int}}, Nmax::Int)

This function is an implementation of a preconditioned CG Method (Polak-Ribière variant of CG), using SAP (pcdSAP function) as a preconditioner 

# Arguments
- **K::SparseMatrixCSC{Float64, Int64} :** K is the sparse stiffness matrix of our inhomogeneous Helmholtz equation
- **F::Matrix{Float64} :** F is the final right hand side of the equation
- **tol::Float64 :** This is the tolerance of the solver (ex. 10^-6) 
- **K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}} :** This contains the stiffness matrices for each subdomain factorized
- **idx::Vector{Vector{Int}} :** The 3 subdomains node vectors with their overlaps
- **Nmax::Int :** The maximum number of iterations

# Returns
 - **q :** The solution of the system K * q = F
"""
function pcdCG(K::SparseMatrixCSC{Float64, Int64}, 
               F::Matrix{Float64}, 
               tol::Float64, 
               K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}},
               idx::Vector{Vector{Int}}, Nmax::Int
               )

    M = size(K, 1);
    q = zeros(M[1]);

    # Initial Residual
    r = F - K * q;
    norm0 = norm(r); # Norm of initial r (r0)

    # Initial preconditioned residual (z0 = M^-1 * r0)
    z = pcdSAP(K, r, K_sub, idx);
    p = z;

    for i in 1:Nmax
        w = K * p;

        a_denom = dot(p, w);
        @assert abs(a_denom) > eps() "pcdCG: Broke down in alpha division (division by 0 = Inf)"
        a = dot(r, z) / a_denom;
        q = q + a*p;
        r_new = r - a*w;

        current_norm = norm(r_new) / norm0;
        println("Iteration $i: ||r|| / ||r0|| = $current_norm");

        if current_norm < tol
            println("pcdCG converged with ", i, " iterations");
            return q
        end

        # New preconditioned residual ( z_k+1 = M^-1 * r_k+1)
        z_new = pcdSAP(K, r_new, K_sub, idx);

        denom2 = dot(r, z);
        @assert abs(denom2) > eps() "pcdCG: Broke down in beta division (division by 0 = Inf)"
        beta = dot(r_new, z_new) / denom2;

        p = z_new + beta*p;

        # Update r and z 
        r = r_new;
        z = z_new;
    end

    println("pcdCG did not converge after $Nmax iterations. Final residual: ", norm(r) / norm0)
    return q
end
