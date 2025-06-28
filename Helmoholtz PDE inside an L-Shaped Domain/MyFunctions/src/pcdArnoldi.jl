"""
    pcdArnoldi(V::AbstractMatrix{Float64}, H::AbstractMatrix{Float64}, 
                effective_operator_apply::Function, m_krylov::Int)

# Arguments
- **V::AbstractMatrix{Float64} :** Orthonormal basis (n x (m_krylov+1))
- **H::AbstractMatrix{Float64} :** Hessenberg matrix (m_krylov+1 x m_krylov)
- **effective_operator_apply::Function :** Function to apply the effective operator
- **m_krylov::Int :** Desired Krylov subspace dimension

# Returns
- **k_actual :** Actual achieved dimension of the Krylov subspace constructed during that particular Arnoldi process
"""
function pcdArnoldi(
    V::AbstractMatrix{Float64},      
    H::AbstractMatrix{Float64},
    effective_operator_apply::Function,
    m_krylov::Int                              
    )
    k_actual = m_krylov;

    for j = 1:m_krylov
        w = effective_operator_apply(V[:, j]); # w = EffectiveOp * v_j

        # Modified Gram-Schmidt
        for i = 1:j
            H[i, j] = dot(V[:, i], w);
            w .-= H[i, j] .* V[:, i];
        end
        
        H[j+1, j] = norm(w);

        if abs(H[j+1, j]) < 1e-12 # Breakdown tolerance
            k_actual = j;
            break
        end
        # Only fill V[:,j+1] if we are not at the last iteration of this m_krylov cycle
        if j <= m_krylov # Check to prevent writing V[:,m_krylov+2] if m_krylov was used in V size
             V[:, j+1] = w / H[j+1, j];
        end
    end
    return k_actual
end