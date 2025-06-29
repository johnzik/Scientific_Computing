# Ihomogenous Helmoholtz PDE Solutions inside an L-Shaped Domain 

- This repository contains the solution of an **Ihomogenous Helmoholtz PDE inside an L-Shaped Domain** written in Julia.


- The PDE is solved using various **numerical methods** as a project for the course "Advanced Scientific Calculation" at my university, Democritus University of Thrace.


- The methods used are:
    1) LU Decomposition
    2) Conjugate Gradient (CG)
    3) Generalized minimal residual method (GMRes)
    4) Domain Decomposition as a preconditioner (pcd)
    5) A pcd CG
    6) A pcd GMRes
    7) Multigrid as a preconditioner to CG
 
## Installation and Run
1) Make sure you have Julia programming language installed.
2) Clone this repository or download: Main.jl, Manifest.toml, Project.toml and the MyFunctions folder and put the all in the same directory.
3) Run Main.jl.  When it runs it will automatically check for the necessary packages, if something isn't installed or is in an unsupported version it will automatically be installed.
4) After a succeful run you should see the method output in your terminal (different methods have different outputs) and 3 plots:
    1) The meshed domain
    2) The solution
    3) The overlapping subdomains (this is used only by pcdCG, pcdGMRes, MGpcdCG) 
