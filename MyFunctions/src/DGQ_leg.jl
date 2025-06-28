"""
    DGQ_leg(f::Function, n::Integer)

# Arguments
- f::Function : The function f(x,y) you want to double integrate from -1 to 1 (This calculates the Legendre GQ)
- n::Integer : The degree of the Gausian Quadrature. With n-point GQ you can solve for polynomials of degree 2n-1 or less

# Returns
- res : The aproximate result of the Double Integration of f(x,y) from -1 to 1
"""
function DGQ_leg(f::Function, n::Integer) # Up to 5 point Double Gaussian Quadrature (Legendre)

    # Gaussian Quadrature points and weights table
    GQ_Table = [
        ([0.0], [2.0]),
        ([+1/sqrt(3), - 1/sqrt(3)], [1.0, 1.0]),
        ([0, +sqrt(3/5), -sqrt(3/5)], [8/9, 5/9, 5/9]),
        ([+sqrt(3/7-(2/7)*sqrt(6/5)), -sqrt(3/7-(2/7)*sqrt(6/5)), +sqrt(3/7+(2/7)*sqrt(6/5)), -sqrt(3/7+(2/7)*sqrt(6/5))]
            ,[(18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36, (18-sqrt(30))/36]),
        ([0, +(1/3)*sqrt(5-2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)), +(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5+2*sqrt(10/7))]
            ,[128/225, (322+13*sqrt(70))/900, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900, (322-13*sqrt(70))/900])
    ]

    # Points and Weights initialization
    p, w = GQ_Table[n]; 

    # Double Gaussian Quadrature calculation
    res = 0.0;
    for i = 1:n
        for j = 1:n
            res += + w[i] * w[j] * f(p[i], p[j]);
        end
    end

    return res
end