include("cg.jl")
include("actionA.jl")
include("actionB.jl")
include("actionBT.jl")
include("data.jl")
include("pcg.jl")
using LinearAlgebra


function inexact_uzawa(n,v1,v2,tol2=1e-6, tol=1e-8, max_iter=10,pcg_max=10,pre_max=2,alpha=1.0,silent=true)
    u = zeros(n-1,n)
    v = zeros(n,n-1)
    p = zeros(n,n)
    f,g = gen_rhs(n)
    res0 = sqrt(norm(f)^2+norm(g)^2 )
    for i = 1:max_iter
        if !silent
            println("Beginning Inexact Uzawa Iteration: $i ...")
        end
        u, v = pcg_inner_loop(p,f,g,tol2,v1,v2,pcg_max,pre_max,silent)
        p = p + alpha*actionBT(u,v)
        bp1, bp2 = actionB(p)
        res1 = actionA_u(u) + bp1 - f
        res2 = actionA_v(v) + bp2 - g
        res3 = actionBT(u,v)
        rescur = sqrt(norm(res1)^2+norm(res2)^2+norm(res3)^2)
        relae = rescur/res0
        if !silent
        println("Inexact Uzawa Iteration: $i, Residual: $relae")
        end
        if relae < tol
            if !silent
            println("Ineaxt Uzawa Converged in ", i, " iterations")
            end
            pde_er = pde_error_sol(u,v)
            return pde_er,i, u, v, p
        end

       
    end
    if !silent
    println("Inexact Uzawa Did not converge after $max_iter iterations.")
    end
    pde_er = pde_error_sol(u,v)
    return pde_er,max_iter, u, v, p
end