"""
    Arnoldi_ModGM(V::AbstractMatrix{Float64}, H::AbstractMatrix{Float64}, A::SparseMatrixCSC, m_krylov::Int)

This function is used by GMRes to compute the Hessenberg matrix H, the Orthonormal basis vectors and the effective Krylov subspace dimension m

# Arguments
- **V::AbstractMatrix{Float64} :** The Orthonormal basis from GMRes function
- **H::AbstractMatrix{Float64} :** The Hessenberg matrix
- **A::SparseMatrixCSC :** This would be our stiffness matrix K
- **m_krylov::Int :** Krylov subspace dimension m

# Returns
- **k_arnoldi :** The new Krylov subspace dimension m
"""
function Arnoldi_ModGM(
    V::AbstractMatrix{Float64}, 
    H::AbstractMatrix{Float64},
    A::SparseMatrixCSC,
    m_krylov::Int)

    k_arnoldi = m_krylov;
    
    for j = 1:m_krylov
    
        w = A * V[:, j];

        # Modified Gram-Schmidt
        for i = 1:j
            H[i, j] = dot(V[:, i], w);
            w .-= H[i, j] .* V[:, i];
        end

        H[j+1, j] = norm(w);

        if abs(H[j+1, j]) < 1e-12 # Breakdown tolerance
            k_arnoldi = j;
            break
        end

        V[:, j+1] = w / H[j+1, j];
    end

    return k_arnoldi
end