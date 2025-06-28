"""
    pcdGMRes(A::SparseMatrixCSC, F_matrix::AbstractMatrix{Float64}, tol::Float64,
            Nmax_restarts::Int, K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}},
            idx::Vector{<:Vector{<:Integer}}; restart_len::Int = min(30, size(A,1)),
            x0_vec::AbstractVector{Float64} = zeros(Float64, size(A,1))
            )

This function an implementation of a preconditioned GMRes method, using SAP (pcdSAP function) as a preconditioner 

# Arguments
- **A::SparseMatrixCSC :** The sparse stiffness matrix of our PDE
- **F_matrix::AbstractMatrix{Float64} :** The RHS of our system
- **tol::Float64 :** This is the tolerance of the solver (ex. 10^-6) 
- **Nmax_restarts::Int :** This is the number of maximum allowed iterations
- **K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}} :** This contains the stiffness matrices for each subdomain factorized
- **idx::Vector{<:Vector{<:Integer}} :** The 3 subdomains node vectors with their overlaps
- **restart_len::Int = min(30, size(A,1)) :** Krylov subspace dimension m
- **x0_vec::AbstractVector{Float64} = zeros(Float64, size(A,1)) :** Initial guess

# Returns
- **x :** The solution of the linear system K * x = F
"""
function pcdGMRes(
    A::SparseMatrixCSC,
    F_matrix::AbstractMatrix{Float64},
    tol::Float64,
    Nmax_restarts::Int,
    K_sub::Vector{SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}},
    idx::Vector{<:Vector{<:Integer}};
    restart_len::Int = min(30, size(A,1)), # Krylov subspace dimension m
    x0_vec::AbstractVector{Float64} = zeros(Float64, size(A,1)) # Initial guess
    )

    n = size(A, 1);
    x = copy(x0_vec); # Initial guess vector
    b_vec = vec(F_matrix); # Ensure RHS is a 1D vector

    total_arnoldi_steps = 0;
    
    initial_r = b_vec - A * x;
    norm0 = norm(initial_r);

    println("Starting GMRES_with_pcdSAP: Max Restarts = $Nmax_restarts, Restart Length = $restart_len, Tolerance = $tol");
    println("Restart 0 (Initial state): RelRes (to initial) approx 1.0, AbsRes = $norm0");

    # Define the preconditioner application function using the wrapper function
    apply_preconditioner = r_v -> pcdSAP_vector_interface_wrapper(A, r_v, K_sub, idx);
    
    # Define the effective operator for Arnoldi: M^-1 * A * v
    effective_A_operator = v_for_mv -> apply_preconditioner(A * v_for_mv);

    for restart_iter = 1:Nmax_restarts
        # Current unpreconditioned residual r_k = b - A*x_k
        upcd_rk = b_vec - A * x;
        
        # Starting vector for Arnoldi: r_k = M^-1 * upcd_rk
        rk = apply_preconditioner(upcd_rk);
        beta = norm(rk);

        if abs(beta) < 1e-14 # Break if preconditioned residual is (near) zero
            println("GMRES_with_pcdSAP: Preconditioned residual norm ($beta) is effectively zero at restart $restart_iter.");
            break 
        end

        Vm = zeros(Float64, n, restart_len + 1); # Initialize Vm
        H = zeros(Float64, restart_len + 1, restart_len); # Hessenberg
        Vm[:, 1] = rk / beta; # Normalized vector

        m_eff = pcdArnoldi(Vm, H, effective_A_operator, restart_len);
        total_arnoldi_steps += m_eff;
        
        if m_eff == 0; break; end 
        
        Hm = H[1:(m_eff+1), 1:m_eff];
        e1 = zeros(Float64, m_eff + 1); e1[1] = 1.0;
        rhs_ls = beta .* e1;
        
        ym = Hm \ rhs_ls; # Least squares solve for ym (length m_eff)
        
        x .+= Vm[:, 1:m_eff] * ym; # Update solution: x_new = x_old + V_m * y_m
        
        # Convergence check -> Relative residual to initial unpreconditioned residual
        r = b_vec - A * x;
        relative_res = norm(r) / norm0;
        
        println("GMRES_with_pcdSAP Restart #$restart_iter (Arnoldi dim $m_eff): RelResToInitial = $relative_res, AbsRes = $(norm(r))");

        if norm(r) < tol || relative_res < tol
            println("GMRES_with_pcdSAP converged with $total_arnoldi_steps iterations");
            break
        end
        if restart_iter == Nmax_restarts
            println("GMRES_with_pcdSAP reached max restarts ($Nmax_restarts).");
        end
    end

    return x
end