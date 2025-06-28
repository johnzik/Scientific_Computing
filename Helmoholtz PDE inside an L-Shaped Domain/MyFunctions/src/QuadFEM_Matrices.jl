"""
    QuadFEM_Matrices(GO::Vector{Matrix{Int64}}, f::Function, k::Float64)

# Arguments
- GO::Vector{Matrix{Int64}} : The Global Ordering matrix which contains the 3 seperate subdomains matrices meshed (use Lmesh to create it)
- f::Function : The right hand side function of the PDE
- k::Float64 : k is the wave number (The bigger it is we expect more peaks to appear)

# Returns
- K : K is the sparse stiffness matrix of our inhomogeneous Helmholtz equation 
- F_global : F_global is the final right hand side of the equation. Now use a solver to solve the system K*u=F_global for u(x,y,t)
"""
function QuadFEM_Matrices(GO::Vector{Matrix{Int64}}, f::Function, k::Float64)

    rowNodes = length(GO[1][1,:]);
    colNodes = length(GO[1][:,1]);
    nop = 3*rowNodes^2 - 2*rowNodes; # Number of points
    rowElem = rowNodes - 1;
    colElem = colNodes - 1;
    noe = rowElem*colElem*3; # Number of elements
    step1 = 1/(rowNodes - 1); # Get step from the GO matrix

    # ---- Create local 2 global mapping ----------------------------------------------------------------------

    l2g = zeros(Int, noe, 4); 

    idx = 1; # Index
    for i = 0:rowElem-1
        for j = 1:colElem
            l2g[idx, :] = [GO[1][end-i, j], GO[1][end-i, j+1], GO[1][end-i-1, j+1], GO[1][end-i-1, j]];
            l2g[idx + Int(noe/3), :] = [GO[2][end-i, j], GO[2][end-i, j+1], GO[2][end-i-1, j+1], GO[2][end-i-1, j]];
            l2g[idx + Int(2*noe/3), :] = [GO[3][end-i, j], GO[3][end-i, j+1], GO[3][end-i-1, j+1], GO[3][end-i-1, j]];
            idx = idx + 1;
        end
    end

    # ---- Point Coords ---------------------------------------------------------------------------------------

    coords = zeros(nop,2);

    # Coordinates of points inside GO[1] & GO[2]
    for i = 1:rowNodes
        xvals = collect(-1:step1:1); # x values from -1 to 1 with step1
        yvals = (-1 + step1*(i - 1))*ones(rowNodes*2-1,1); # rowNodes*2-1 cause there is 1 common point, increase from -1 
                                                           # with step1 once every iteration

        start_idx = (i - 1)*(rowNodes*2-1) + 1;
        end_idx = (rowNodes*2-1)*i;
        coords[start_idx:end_idx, :] = hcat(xvals, yvals);
    end

    # Coordinates of points inside GO[3]
    for i = 2:rowNodes
        xvals = collect(-1:step1:0); # x values from -1 to 1 with step1
        yvals = (0 + step1*(i - 1))*ones(rowNodes,1); # rowNodes*2-1 cause there is 1 common point, increase from -1 
                                                      # with step1 once every iteration

        start_idx = GO[3][end-i+1,1];
        end_idx = GO[3][end-i+1,end];
        coords[start_idx:end_idx, :] = hcat(xvals, yvals);
    end

    # ---- Set Variables & Iterate and Populate ---------------------------------------------------------------

    # Interpolation functions
    N = [
        (ksi, h) -> (1/4)*(1-ksi)*(1-h);
        (ksi, h) -> (1/4)*(1+ksi)*(1-h);
        (ksi, h) -> (1/4)*(1+ksi)*(1+h);
        (ksi, h) -> (1/4)*(1-ksi)*(1+h); 
    ]

    # Assemble Stiffness Matrix & Right Hand Side 
    Me = zeros(4, 4); # Element Me
    Te = zeros(4, 4); # Element Te

    ia = zeros(Int, 16*noe); # Row Sparse idx
    ja = zeros(Int, 16*noe); # Column Sparse idx
    va = zeros(Float64, 16*noe); # Value Sparse idx

    fe = zeros(4); # Element f
    F_global = zeros(nop, 1); # Global F

    c = 1; # idx for Sparse construction
    for e = 1:noe
        # Point Coordinates of each element
        xe = coords[Int.(l2g[e,:]), 1]; 
        ye = coords[Int.(l2g[e,:]), 2];

        # Jacobian matrix values quadrilateral elements 
        J11 = (ksi, h) -> (1/4)*(-(1-h)*xe[1] + (1-h)*xe[2] + (1+h)*xe[3] - (1+h)*xe[4]);
        J12 = (ksi, h) -> (1/4)*(-(1-h)*ye[1] + (1-h)*ye[2] + (1+h)*ye[3] - (1+h)*ye[4]);
        J21 = (ksi, h) -> (1/4)*(-(1-ksi)*xe[1] - (1+ksi)*xe[2] + (1+ksi)*xe[3] + (1-ksi)*xe[4]);
        J22 = (ksi, h) -> (1/4)*(-(1-ksi)*ye[1] - (1+ksi)*ye[2] + (1+ksi)*ye[3] + (1-ksi)*ye[4]);
        detJ = (ksi, h) -> J11(ksi,h)*J22(ksi,h) - J12(ksi,h)*J21(ksi,h); # Determinant of J matrix

        # Interpolation functions derivatives
        dNx = [
            (ksi, h) -> (1/(4*detJ(ksi, h))) * (-J22(ksi, h)*(1-h) + J12(ksi, h)*(1-ksi));
            (ksi, h) -> (1/(4*detJ(ksi, h))) * ( J22(ksi, h)*(1-h) + J12(ksi, h)*(1+ksi));
            (ksi, h) -> (1/(4*detJ(ksi, h))) * ( J22(ksi, h)*(1+h) - J12(ksi, h)*(1+ksi));
            (ksi, h) -> (1/(4*detJ(ksi, h))) * (-J22(ksi, h)*(1+h) - J12(ksi, h)*(1-ksi));
        ]

        dNy = [
            (ksi, h) -> (1/(4*detJ(ksi, h))) * ( J21(ksi, h)*(1-h) - J11(ksi, h)*(1-ksi));
            (ksi, h) -> (1/(4*detJ(ksi, h))) * (-J21(ksi, h)*(1-h) - J11(ksi, h)*(1+ksi));
            (ksi, h) -> (1/(4*detJ(ksi, h))) * (-J21(ksi, h)*(1+h) + J11(ksi, h)*(1+ksi));
            (ksi, h) -> (1/(4*detJ(ksi, h))) * ( J21(ksi, h)*(1+h) + J11(ksi, h)*(1-ksi));
        ]

        for i = 1:4 # 1:4 because we have 4 nodes per element
            for j = 1:4 
                # Assemble Stiffness Matrix
                Mfunc = (ksi, h) -> (dNx[i](ksi, h)*dNx[j](ksi, h) + dNy[i](ksi, h)*dNy[j](ksi, h)) * detJ(ksi, h);
                Me[i,j] = - DGQ_leg(Mfunc, 3);

                Tfunc = (ksi, h) -> k^2 * N[i](ksi, h) * N[j](ksi, h) * detJ(ksi, h);
                Te[i,j] = DGQ_leg(Tfunc, 3);
            end

            # Assemble right hand side
            xkh = (ksi, h) -> xe[1]*N[1](ksi, h) + xe[2]*N[2](ksi, h) + xe[3]*N[3](ksi, h) + xe[4]*N[4](ksi, h);
            ykh = (ksi, h) -> ye[1]*N[1](ksi, h) + ye[2]*N[2](ksi, h) + ye[3]*N[3](ksi, h) + ye[4]*N[4](ksi, h);
            ffunc = (ksi, h) -> N[i](ksi, h)*(f(xkh(ksi,h), ykh(ksi,h))) * detJ(ksi, h);
            fe[i] = DGQ_leg(ffunc ,3);

            # Global f
            F_global[Int.(l2g[e, i])] = F_global[Int.(l2g[e, i])] + fe[i];
        end

        # Global Sparse Striffness Matrix Construction
        for i = 1:4
            for j = 1:4
                ia[c] = l2g[e,i];
                ja[c] = l2g[e,j];
                va[c] = Me[i,j] + Te[i,j];
                c += 1;
            end
        end

    end

    K = sparse(ia, ja, va, nop, nop);

    # ---- Find & Enforce Dirichlet Boundary Conditions -------------------------------------------------------
   
    # Find the nodes
    boundary_nodes = findall(
        (coords[:,1] .== -1) .| 
        (coords[:,1] .== 1) .| 
        (coords[:,2] .== -1) .| 
        (coords[:,2] .== 1) .| 
        ((coords[:,1] .== 0) .& (coords[:,2] .>= 0)) .| 
        ((coords[:,2] .== 0) .& (coords[:,1] .>= 0))
    );

    # Enforce the Boundary Conditions
    for i in boundary_nodes
        cols, _ = findnz(K[i, :])  # Ignore 2nd output of findnz (which is the actual non-zero values themselves)
                                   # and keep only the columns indices of the non-zero elements in row i
        for j in cols
            K[i, j] = 0.0
        end
        K[i, i] = 1.0
        F_global[i] = 0.0
    end

    return K, F_global, coords;
end