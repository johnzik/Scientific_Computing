struct MGLevel
    A::SparseMatrixCSC{Float64, Int64} # Stiffness matrix
    R::Union{SparseMatrixCSC{Float64, Int64}, Nothing} # Restriction (nothing for coarsest)
    P::Union{SparseMatrixCSC{Float64, Int64}, Nothing} # Prolongation (nothing for coarsest)
    nodes_original_indices::Union{Vector{Int}, Nothing}
    coords::Union{Matrix{Float64}, Nothing}             
end

struct MGHier
    levels::Vector{MGLevel}
end