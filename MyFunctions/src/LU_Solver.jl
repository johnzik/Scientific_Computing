"""
    LU_Solver(K::SparseMatrixCSC{Float64, Int64}, F_global::Matrix{Float64})

This function solves a Sparse Linear System, K * x = F, using LU Factorization and outputs
the solution   

# Arguments
- **K::SparseMatrixCSC{Float64, Int64} :** The Sparse Global Stiffness Matrix of our PDE
- **F_global::Matrix{Float64} :** The Global RHS of our PDE

# Returns
- **X_solution :** This is result after solving the system
"""
function LU_Solver(K::SparseMatrixCSC{Float64, Int64}, F_global::Matrix{Float64})
    F = lu(K);
    L = F.L; # Unit Lower Triangular Sparse Matrix
    U = F.U; # Upper Triangular Sparse Matrix
    p = F.p; # Row permutation vector
    q = F.q; # Column permutation vector
    Rs = F.Rs; # Row scaling vector

    n = size(K, 1);
    num_rhs = size(F_global, 2);
    X_solution = Matrix{Float64}(undef, n, num_rhs);

    # This single working vector will be transformed in-place:
    # F_col -> b' -> y_col -> z_col
    work_vector = Vector{Float64}(undef, n);

    for k_rhs = 1:num_rhs
        # Step 1) Prepare RHS: work_vector = b' = Rs * P * b
        for i = 1:n
            work_vector[i] = Rs[i] * F_global[p[i], k_rhs];
        end

        # Step 2) We have L * U * z = b' = work_vector
        # Let U * z = y -> Solve L * y = b' for y (Forward Substitution)
        for j = 1:n  # Iterate through columns of L
            val_y_j = work_vector[j]; # This y_j is final (as L_jj=1)
            if val_y_j != 0.0
                for k_idx = L.colptr[j]:(L.colptr[j+1]-1)
                    i_row = L.rowval[k_idx];
                    if i_row > j # Element L_ij is below the diagonal
                        work_vector[i_row] -= L.nzval[k_idx] * val_y_j;
                    end
                end
            end
        end
        # Now, work_vector = y

        # Step 3) U * z = y = work_vector -> Solve for z (Backward Substitution)
        for j = n:-1:1  # Iterate backwards through columns of U
            
            # Find U_jj (diagonal element of U in column j)
            U_jj = 0.0;

            for k_idx = U.colptr[j]:(U.colptr[j+1]-1)
                if U.rowval[k_idx] == j
                    U_jj = U.nzval[k_idx];
                    break
                end
            end

            current_sum_val = work_vector[j]; # This is y_j - sum(U_jk * z_k for k>j)
            
            if U_jj == 0.0
                work_vector[j] = NaN; # Matrix is singular
            else
                work_vector[j] = current_sum_val / U_jj;
            end
            
            val_z_j = work_vector[j]; # This is the final value for z_j

            if val_z_j != 0.0 && !isnan(val_z_j)
                for k_idx = U.colptr[j]:(U.colptr[j+1]-1)
                    i_row = U.rowval[k_idx];
                    if i_row < j # Element U_ij is above the diagonal
                        work_vector[i_row] -= U.nzval[k_idx] * val_z_j;
                    end
                end
            end
        end
        # Now, work_vector = z

        # Step 4) x = Q * z -> Apply inverse column permutation -> Get solution x
        current_X_col = view(X_solution, :, k_rhs)
        for i = 1:n
            current_X_col[q[i]] = work_vector[i];
        end
    end

    return X_solution
end