include("actionA.jl")
include("vcycle.jl")
include("actionBT.jl")
include("actionB.jl")


function pcg_inner_loop(P,f,g,tol,v1,v2,pcg_max_iter=10,pre_max=2,silent=false)
    bp1, bp2 = actionB(P)
    rhs1 = f - bp1
    rhs2 = g - bp2
    xu = zeros(size(rhs1))
    pu = zeros(size(rhs1))
    ru = zeros(size(rhs1))
    zu = zeros(size(rhs1))

    xu2 = zeros(size(rhs1))
    pu2 = zeros(size(rhs1))
    ru2 = zeros(size(rhs1))
    zu2 = zeros(size(rhs1))


    xu, pu, ru, zu = pcg_init_u(rhs1, xu, v1, v2, pre_max, silent)

    xv = zeros(size(rhs2))
    pv = zeros(size(rhs2))
    rv = zeros(size(rhs2))
    zv = zeros(size(rhs2))

    xv2 = zeros(size(rhs2))
    pv2 = zeros(size(rhs2))
    rv2 = zeros(size(rhs2))
    zv2 = zeros(size(rhs2))

    xv, pv, rv, zv = pcg_init_v(rhs2, xv, v1, v2, pre_max, silent)

    for iter in 1:pcg_max_iter
        #check convergence 
        expr1 = actionA_u(xu) - rhs1
        expr2 = actionA_v(xv) - rhs2
        expr3 = actionBT(xu, xv)
        criterion = sqrt(norm(expr1)^2+norm(expr2)^2)/norm(expr3)
        if !silent
            println("PCG inner loop Iteration $iter: Relative residual = $criterion")
        end

        if criterion < tol
             if !silent
                println("PCG inner loop converged after $iter iterations.")
             end 
            return xu, xv
        end

        xu2, pu2, ru2, zu2 = pcg_update_u(xu, pu, ru, zu, v1, v2, pre_max, silent)
        xv2, pv2, rv2, zv2 = pcg_update_v(xv, pv, rv, zv, v1, v2, pre_max, silent)

        xu = xu2 
        pu = pu2
        ru = ru2
        zu = zu2
        xv = xv2
        pv = pv2
        rv = rv2
        zv = zv2
    end
    if !silent
        println("PCG inner loop did not converge after $pcg_max_iter iterations.")
    end
    return xu, xv
    



end



function pcg_update_u(x,p,r,z,v1,v2,pre_max=2,silent=false)
    xnew = zeros(size(x))
    rnew = zeros(size(r))
    znew = zeros(size(z))
    pnew = zeros(size(p))
    rsold = sum(r .* z)
    AP = actionA_u(p)
    alpha = rsold / sum(p .* AP)
    xnew = x + alpha .* p
    rnew = r - alpha .* AP
    znew, it = precon_Au(rnew, v1, v2, pre_max)
    
    rsnew = sum(rnew .* znew)
    beta = rsnew / rsold
    pnew = znew + beta .* p
    return xnew, pnew, rnew, znew

end

function pcg_update_v(x,p,r,z,v1,v2,pre_max=2,silent=false)
    xnew = zeros(size(x))
    rnew = zeros(size(r))
    znew = zeros(size(z))
    pnew = zeros(size(p))
    rsold = sum(r .* z)
    AP = actionA_v(p)
    alpha = rsold / sum(p .* AP)
    xnew = x + alpha .* p
    rnew = r - alpha .* AP
    znew, it = precon_Av(rnew, v1, v2, pre_max)
   
    rsnew = sum(rnew .* znew)
    beta = rsnew / rsold
    pnew = znew + beta .* p
    return xnew, pnew, rnew, znew

end


function pcg_init_u(rhs, x0, v1, v2,pre_max=3,silent = false)
    #return x ,p, r, z
    r = rhs - actionA_u(x0)
    z = zeros(size(r))
    z, it = precon_Au(r, v1, v2, pre_max)
    
    p = copy(z)
    return x0, p, r, z
end


function pcg_init_v(rhs, x0, v1, v2,pre_max=3,silent = false)
    #return x ,p, r, z
    r = rhs - actionA_v(x0)
    z = zeros(size(r))
    z, it = precon_Av(r, v1, v2, pre_max)
    
    p = copy(z)
    return x0, p, r, z
end







function pcg_solver_u(R, u0, v1, v2,pre_max=2, maxiter=100, tol=1e-6, tol2=1e-3, silent=false)
    # Initial guess for U is provided by u0
    U = u0
    # Calculate initial residual
    Residual = R - actionA_u(U)
    # Apply preconditioner to the initial residual
    Z, it = precon_Au(Residual, tol2, v1, v2, pre_max, true)
    if !silent
        println("vcycle iterations: ", it)
    end
    P = copy(Z)
    rsold = sum(Residual .* Z)
    success = false
    for iter in 1:maxiter
        AP = actionA_u(P)
        alpha = rsold / sum(P .* AP)
        U .+= alpha .* P
        Residual .-= alpha .* AP
        # Apply preconditioner to the updated residual
        Z, it = precon_Au(Residual, tol2, v1, v2, pre_max, true)
        if !silent
            println("vcycle iterations: ", it)
        end

        rsnew = sum(Residual .* Z)
        
        # Check convergence
        ind = norm(actionA_u(U) - R)/norm(R)
        if !silent
            println("PCG-U Iteration $iter: Relative residual = $ind")
        end

        if ind < tol
            if !silent
                println("PCG-U Convergence achieved after $iter iterations.")
            end
            success = true
            break
        end

        beta = rsnew / rsold
        P = Z .+ beta .* P
        rsold = rsnew
    end
    
    if !success && !silent
        println("PCG - U Maximum iterations reached without convergence.")
    end

    return U
end



