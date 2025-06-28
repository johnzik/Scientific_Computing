function ndgrid(x::AbstractVector, y::AbstractVector)
    X = repeat(x', length(y), 1)  # Transpose x, repeat across rows
    Y = repeat(y, 1, length(x))   # Repeat y across columns
    return X, Y
end