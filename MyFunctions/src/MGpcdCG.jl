function MGpcdCG(
    K::SparseMatrixCSC, 
    F::AbstractMatrix{Float64},
    tol::Float64, 
    idx::Vector{Vector{Int64}},
    mg_hier::Vector{MGHier},
    Nmax::Int,
    omega::Float64 = 4/5, # Default smoother omega
    nPre::Int = 2, # Default number of pre-smoothing iterations
    nPost::Int = 1 # Default number of post-smoothing iterations
    )
    F_vec = vec(F);

    nop = size(K, 1);
    q = zeros(Float64, nop);

    r = F_vec - K * q; # Initial residual (vector)
    norm0 = norm(r);
    if norm0 < 1e-14
        println("MGpcdCG: Initial residual is zero. Solution is the initial guess.")
        return q
    end
    
    # Initial preconditioned residual (z0 = M^-1 * r0)
    z = MGpcdSAP(K, r, idx, mg_hier, nPre, nPost, omega);
    p = copy(z);

    rz_old = dot(r, z) # For Polak-Ribière beta

    println("Starting MGpcdCG: Max Iterations = $Nmax, Tolerance = $tol")
    println("Iteration 0: ||r||/||r0|| = 1.0")
    
    for iter in 1:Nmax
        w = K * p;
        
        a_denom = dot(p, w);
        @assert abs(a_denom) > eps() "MGpcdCG: Broke down in alpha division (division by 0 = Inf)"
        a = rz_old / a_denom;
        q .+= a .* p;
        r .-= a .* w; # r is now r_{k+1}
        
        current_rel_res = norm(r) / norm0
        println("Iteration $iter: ||r|| / ||r0|| = $current_rel_res")

        if current_rel_res < tol
            println("MGpcdCG converged in $iter iterations.");
            return q
        end
        
        # New preconditioned residual (z_{k+1} = M^-1 * r_{k+1})
        z_new = MGpcdSAP(K, r, idx, mg_hier, nPre, nPost, omega);
        
        @assert abs(rz_old) > eps() "MGpcdCG: Broke down in beta division (division by 0 = Inf)"
        rz_new = dot(r, z_new);
        beta = rz_new / rz_old; # Polak-Ribière: (r_{k+1}^T z_{k+1}) / (r_k^T z_k)
        
        p = z_new .+ beta .* p;

        rz_old = rz_new;
    end
    
    println("MGpcdCG did not converge after $Nmax iterations. Final relative residual: ", norm(r) / norm0)
    return q
end