"""
    Conj_Grad(A, B, tol)

Return the solution of 'A * x = B' using the Conjugate Gradient Method (CG)
"""
function Conj_Grad(A, B, tol, Nmax)
    M = size(A)
    q = zeros(M[1])

    # Residual
    r = B - A * q
    p = r

    # Norm of rbs
    norm0 = norm(r)

    if norm(r) > tol * norm0
        for i in 1:Nmax
            denom = dot(p, A*p);
            @assert abs(denom) > eps() "CG: Broke down in alpha division (division by 0 = Inf)"
            a = dot(r, r) / denom; # dot(a, b) computes (a' * b)
            q += a * p;

            temp_r = r;
            r = r - a * A * p;

            if norm(r) < tol * norm0
                println("CG converged with ", i, " iterations");
                return q
            end

            denom2 = dot(temp_r, temp_r);
            @assert abs(denom2) > eps() "CG: Broke down in beta division (division by 0 = Inf)"
            beta = dot(r, r) / denom2;
            p = r + beta * p;
        end
    end

    println("CG did not converge after $Nmax iterations. Final residual: ", norm(r) / norm0)
    return q
end