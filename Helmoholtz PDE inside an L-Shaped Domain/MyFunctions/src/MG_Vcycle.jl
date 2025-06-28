function MG_Vcycle(
    b::Vector{Float64},
    mg::MGHier,
    nPre::Int,
    nPost::Int,
    omega::Float64,
    x0_guess::Vector{Float64}=zeros(Float64, length(b)),
    lvl_idx::Int = 1,
    )

    lev_data = mg.levels[lvl_idx];
    A = lev_data.A;

    # Coarsest level check
    if lvl_idx == length(mg.levels) || lev_data.R == nothing
        if size(A,1) > 0
            return A \ b
        else # A is empty, b should also be empty
            return zeros(Float64,0)
        end
    end

    # Initialize Solution
    x = copy(x0_guess);

    # Pre-smoothing (Damped Jacobi)
    Dinv = 1.0 ./ diag(A);
    for _ = 1:nPre
        res_smooth_pre  = b - A * x;
        x .+= omega .* (res_smooth_pre .* Dinv);
    end
    # gauss_seidel(x, A, b, nPre, omega_sor=omega);

    # Restrict
    rFine = b - A * x;
    rCoarse = lev_data.R * rFine;

    # Recursive correction
    e_coarse_initial = zeros(Float64, length(rCoarse));
    e_coarse = MG_Vcycle(
        rCoarse, 
        mg,
        nPre,
        nPost,
        omega,
        e_coarse_initial,
        lvl_idx + 1
    );

    # Prolongate and correct
    if !isempty(e_coarse) # Check if coarse error is not empty
        x .+= lev_data.P * e_coarse
    end

    # Post-smoothing
    for _ = 1:nPost
        res_smooth_post = b - A * x;
        x .+= omega .* (res_smooth_post .* Dinv)
    end
    # gauss_seidel(x, A, b, nPost, omega_sor=omega);

    return x
end