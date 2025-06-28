using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate() 
using MyFunctions, BenchmarkTools, SparseArrays, LinearAlgebra;
import PrettyTables;

step1 = 0.01;
GO = Lmesh(step1);

# Define the function as an anonymous function
s = 0.05;
k = 8.0; # If k = 0 then we have the Poisson Equation (Because we have division by zero errors use eps(Float64) instead of 0). The bigger k is we have more peaks
f = (x, y) -> (-1/(2*pi*s^2)) * exp(-((x+0.5).^2 + (y-0.5).^2) / (2*s^2)) + (-1/(2*pi*s^2)) * exp(-((x-0.5).^2 + (y+0.5).^2) / (2*s^2));

println("Building global FEM matrices (k=$k)...");
@time "Global FEM Assembly" K1, F_global1, all_coords = QuadFEM_Matrices(GO, f, k);     

H = 40; # Overlap
println("\nPerforming domain decomposition with overlap H = $H...")
@time "Subdomain Overlap" K_sub, F_sub, idx = SubdomainOverlap_Matrices(GO, K1, F_global1, H);
println("Number of subdomains created: ", length(K_sub));
K_sub_factors = [factorize(K_s) for K_s in K_sub]; # Pre-factorize K_sub to improve performance

# --- Multigrid Setup Parameters ---
num_mg_levels = 2; # Number of levels in hierarchy
min_coarse_size_param = 10; # Min nodes on coarsest grid for hierarchy building
# V-cycle parameters
nPre= 2; # Pre-smoothing iterations
nPost = 2; # Post-smoothing iterations
omega = 4.0/5.0; # Damping for Jacobi (0.9 got the best results)

# Setup the MG hierarchies
mg_hierarchies_list = MGHier_Setup(K_sub, idx, all_coords, num_mg_levels, min_coarse_size_param);

# --- Solving the System ---
tol = 10^-4;
Nmax = 500; # Max iterations
krylov_dim_gmres = 30;

println("\n--- Solving with k=$k ---")

@time "Solving time" begin
    # u = K1\F_global1;
    # u = LU_Solver(K1, F_global1);
    # u = Conj_Grad(K1, F_global1, tol, 1000);
    # u = GMRes(K1, F_global1, tol, max_restarts = Nmax);
    # u = pcdGMRes(K1, F_global1, tol, Nmax, K_sub_factors, idx; restart_len = krylov_dim_gmres);
    # u = pcdCG(K1, F_global1, tol, K_sub_factors, idx, 500);
    u = MGpcdCG(K1, F_global1, tol, idx, mg_hierarchies_list, Nmax, omega, nPre, nPost);
end


# ----------------------------------------------------------------------------
#       Plotting the result
# ----------------------------------------------------------------------------

PlotDomain(all_coords);
PlotDomainOverlaps(all_coords, idx);
PlotSolutionAnimation(GO, u, k);

nothing

