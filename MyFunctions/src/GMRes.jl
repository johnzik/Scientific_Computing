"""
    GMRes(A::SparseMatrixCSC, b_matrix::Matrix{Float64}, tol::Float64; max_restarts::Int = 300,restart_length::Int = min(30, size(A,1)))

Return the solution of 'A * x = B' using the Generalized Minimal Residual method (GMRes) 
"""
function GMRes(
    A::SparseMatrixCSC,
    b_matrix::Matrix{Float64},
    tol::Float64;
    max_restarts::Int = 300,
    restart_length::Int = min(30, size(A,1))
    )
    
    b = vec(b_matrix);
    n = size(A, 1);
    total_iterations = 0;

    x = zeros(Float64, size(A,1));

    # Initial residual
    r = b - A * x;
    norm0 = norm(r);
    
    if norm(r) < tol
        println("GMRES: Initial guess is already a solution (residual norm < tol).");
        return 0
    end
    
    for restart_iter = 1:max_restarts
        
        rk = b - A * x;
        beta = norm(rk);

        V = zeros(Float64, n, restart_length + 1);
        H = zeros(Float64, restart_length + 1, restart_length);
        V[:, 1] = rk / beta;

        m = Arnoldi_ModGM(V, H, A, restart_length);
        total_iterations += m;

        Hm = H[1:(m+1), 1:m];

        e1_rhs = zeros(Float64, m + 1);
        e1_rhs[1] = 1.0;
        rhs_least_squares = beta .* e1_rhs;
        ym = Hm \ rhs_least_squares;

        x .+= V[:, 1:m] * ym;

        r = b - A * x;

        if norm(r) <= tol * norm0 || norm(r) < tol
            println("GMRES converged with $total_iterations iterations");
            break
        end

        if restart_iter == max_restarts
            println("GMRES reached max restarts ($max_restarts).");
        end
    end

    return x
end