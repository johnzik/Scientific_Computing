module MyFunctions
    using GLMakie, SparseArrays, LinearAlgebra;
    # import PrettyTables;

    # --- Struct ---
    include("MGLevel.jl");

    # --- FEM & Domain Decomposition ---
    include("Conj_Grad.jl");
    include("DGQ_leg.jl");
    include("Lmesh.jl");
    include("ndgrid.jl");
    include("pcdCG.jl");
    include("pcdSAP.jl");
    include("QuadFEM_Matrices.jl");
    include("LU_Solver.jl");
    include("SubdomainOverlap_Matrices.jl");
    include("Arnoldi_ModGM.jl");
    include("GMRes.jl");
    include("pcdGMRes.jl");
    include("pcdArnoldi.jl");
    include("pcdSAP_vector_interface_wrapper.jl");

    # --- Multigrid Components ---
    include("MG_Vcycle.jl");
    include("MGpcdSAP.jl");
    include("MGpcdCG.jl");
    include("build_mg_hierarchy.jl");
    include("bilinearRnP.jl");
    include("MGHier_Setup.jl");

    # --- Plotting ---
    include("PlotDomain.jl")
    include("PlotDomainOverlaps.jl");
    include("PlotSolutionAnimation.jl");

    export Conj_Grad, DGQ_leg, Lmesh, ndgrid, pcdCG, pcdSAP, QuadFEM_Matrices, LU_Solver, SubdomainOverlap_Matrices, MGLevel;
    export MG_Vcycle, MGpcdSAP, MGpcdCG, build_mg_hierarchy, MGHier, PlotDomainOverlaps, PlotSolutionAnimation, PlotDomain;
    export Arnoldi_ModGM, GMRes, pcdArnoldi, pcdGMRes, pcdSAP_vector_interface_wrapper, bilinearRnP, MGHier_Setup;
end
