include("cg.jl")
include("actionA.jl")
include("actionB.jl")
include("actionBT.jl")
include("data.jl")
using LinearAlgebra
function uzawa_update(u,v,p,f,g,alpha=1.0)
    n = size(p, 1)
    pnew = zeros(n,n)
    bp1 , bp2 = actionB(p)
    r1 = f - bp1
    r2 = g - bp2
    unew = cg_solver_u(r1,zeros(n-1,n))
    vnew = cg_solver_v(r2,zeros(n,n-1))  
  

    pnew = p + alpha*actionBT(unew,vnew)

    return unew, vnew, pnew

end

function uzawa(n, tol=1e-8, max_iter=10,silent=true)
    u = zeros(n-1,n)
    v = zeros(n,n-1)
    p = zeros(n,n)
    f,g = gen_rhs(n)
    res0 = sqrt(norm(f)^2+norm(g)^2 )
    for i = 1:max_iter
        unew, vnew, pnew = uzawa_update(u,v,p,f,g)
        bp1, bp2 = actionB(pnew)
        res1 = actionA_u(unew) + bp1 - f
        res2 = actionA_v(vnew) + bp2 - g
        res3 = actionBT(unew,vnew)
        rescur = sqrt(norm(res1)^2+norm(res2)^2+norm(res3)^2)
        relae = rescur/res0
        if !silent
        println("Iteration: ", i, " Residual: ", relae)
        end
        if relae < tol
            if !silent
            println("Converged in ", i, " iterations")
            end
            pde_er = pde_error_sol(unew,vnew)
            return pde_er,i, unew, vnew, pnew
        end
        u = unew
        v = vnew
        p = pnew
    end
    if !silent
    println("Uzawa Did not converge")
    end
    pde_er = pde_error_sol(u,v)
    return pde_er,max_iter, u, v, p
end