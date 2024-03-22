function cg_solver_u( R, u0, maxiter=2000, tol=1e-10,silent=false)
    # Initial guess for U is a zero matrix with the same dimensions as R
    U = u0
    # Calculate initial residual
    Residual = R - actionA_u(U)
    P = copy(Residual)
    rsold = sum(Residual .^ 2)
    success = false
    for iter in 1:maxiter
        AP = actionA_u(P)
        alpha = rsold / sum(P .* AP)
        U .+= alpha .* P
        Residual .-= alpha .* AP
        rsnew = sum(Residual .^ 2)
        
        # Check convergence
        if sqrt(rsnew) / norm(R) < tol
            if !silent
            println("CG-U Convergence achieved after $iter iterations.")
            end
            success = true
            break
        end

        beta = rsnew / rsold
        P = Residual .+ beta .* P
        rsold = rsnew
    end
    
    if !success && !silent
    println("CG - U Maximum iterations reached without convergence.")
    end

    return U
end

function cg_solver_v( R, v0, maxiter=2000, tol=1e-10, silent=false)
    # Initial guess for U is a zero matrix with the same dimensions as R
    U = v0
    # Calculate initial residual
    Residual = R - actionA_v(U)
    P = copy(Residual)
    rsold = sum(Residual .^ 2)
    success = false

    for iter in 1:maxiter
        AP = actionA_v(P)
        alpha = rsold / sum(P .* AP)
        U .+= alpha .* P
        Residual .-= alpha .* AP
        rsnew = sum(Residual .^ 2)
        
        # Check convergence
        if sqrt(rsnew) / norm(R) < tol
            if !silent
            println("CG-V Convergence achieved after $iter iterations.")
            end
            success = true
            break
        end

        beta = rsnew / rsold
        P = Residual .+ beta .* P
        rsold = rsnew
    end
    
    if !success && !silent
    println("CG - V Maximum iterations reached without convergence.")
    end

    return U
end
