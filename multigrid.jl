include("actionA.jl")
include("actionB.jl")
include("actionBT.jl")
include("gs.jl")
include("dgs.jl")
include("restrict.jl")
include("lift.jl")
include("cg.jl")
include("basesolution.jl")
include("data.jl")

function vcycle_stokes(N,v1,v2,tol=1e-8,maxiter=10, quiet = false)
    f, g = gen_rhs(N)
    d = zeros(N,N)
    u = zeros(N-1,N)
    v = zeros(N,N-1)
    p = zeros(N,N)
    r1 = zeros(N-1,N)
    r2 = zeros(N,N-1)
    r3 = zeros(N,N)
    u4 = zeros(N-1,N)
    v4 = zeros(N,N-1)
    p4 = zeros(N,N)
    for iter in 1:maxiter
        
       
        u4,v4,p4 = vcycle_dgs(u,v,p,f,g,d,N,v1,v2, true)
        r1,r2,r3, res = cal_res(u4,v4,p4,f,g,d)
        if !quiet
            println("V-cycle iteration $iter, residual: $res")
        end
        if res < tol
            pde_er = pde_error_sol(u4,v4)
            return pde_er,u4,v4,p4, iter
        end
        u = u4
        v = v4
        p = p4
    end
    if !quiet
        println("V-cycle did not converge after $maxiter iterations, residual: $res")
    end
    pde_er = pde_error_sol(u4,v4)
    return pde_er,u4,v4,p4, maxiter



end





function vcycle_dgs(u,v,p, f,g,d, level, v1, v2, quiet = true)
    r1 = zeros(level-1,level)
    r2 = zeros(level,level-1)
    r3 = zeros(level,level)
    err =0.0

    if level == 4
        # Solve exactly at the finest level using ground_truth_4
        if !quiet
            println("Reaching level 4...")
        end
        x1,x2,x3,it = basesolver(f,g,d)
        r1, r2, r3, res = cal_res(x1,x2,x3,f,g,d)
        if !quiet
            println("Level 4, current residual: $res, inner iterations: $it")

        end
        return x1,x2,x3
    else
        # Pre-smoothing with Gauss-Seidel
        if !quiet
            println("Presmoothing level $level ...")
        end

        for _ in 1:v1
            dgs!(u,v,p,f,g,d)
        end

        # Compute residual
        r1,r2,r3, err = cal_res(u,v,p,f,g,d)
        if !quiet
        println("At level $level, After presmoothing,  residual: , $err")
        println("Restricting level $level ...")
        end
        # Restrict residual to coarser grid
        r1coarse = restrict_u_loop(r1)
        r2coarse = restrict_v_loop(r2)
        r3coarse = restrict_p_loop(r3)

        e1coarse = zeros(size(r1coarse))
        e2coarse = zeros(size(r2coarse))
        e3coarse = zeros(size(r3coarse))

        # Recursive call on the coarser grid
        e1solved, e2solved, e3solved = vcycle_dgs(e1coarse,e2coarse,e3coarse, r1coarse,r2coarse,r3coarse, Int(level/2), v1, v2, quiet)

        # Lift the approximation of the error back to the finer grid
        if !quiet
        println("Lifting level $level ...")
        end
        e1fine = lift_u_loop(e1solved)
        e2fine = lift_v_loop(e2solved)
        e3fine = lift_p_loop(e3solved)

        # Correct the approximation on the current grid
        u += e1fine
        v += e2fine
        p += e3fine

        r1,r2,r3, err = cal_res(u,v,p,f,g,d)

        
        if !quiet
        println("At level $level, After eliminating low frequency error, current residual: ", err)
        println("Postsmoothing level $level ...")
        end
        # Post-smoothing with Gauss-Seidel
        for _ in 1:v2
            dgs!(u,v,p,f,g,d)
        end

        r1, r2, r3, resf = cal_res(u,v,p,f,g,d)
        if !quiet
        println("At level $level, After postsmoothing, current residual: ", resf)
        end
        return u,v,p
        
    end
end

function control_update(u,v,p, f,g,d, level, v1, v2, quiet = true)
    r1 = zeros(level-1,level)
    r2 = zeros(level,level-1)
    r3 = zeros(level,level)
    err =0.0
    
    for _ in 1:v1
        dgs!(u,v,p,f,g,d)
    end
    
    for _ in 1:v2
        dgs!(u,v,p,f,g,d)
    end

    r1, r2, r3, res = cal_res(u,v,p,f,g,d)
    return u,v,p, res

end