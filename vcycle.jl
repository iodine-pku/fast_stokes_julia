include("actionA.jl")
include("symmetric_gs.jl")
include("restrict.jl")
include("lift.jl")
include("cg.jl")

function vcycle_Au(u, b, level, v1, v2, quiet = true)
    if level == 4
        # Solve exactly at the finest level using ground_truth_4
        if !quiet
            println("Reaching level 4...")
        end
        sol = cg_solver_u(b, zeros(3,4), 300, 1e-10, quiet)
        res = norm(b - actionA_u(sol))/norm(b)
        if !quiet
            println("Level 4, current residual: ", res)

        end
        return sol
    else
        # Pre-smoothing with Gauss-Seidel
        if !quiet
            println("Presmoothing level $level ...")
        end

        for _ in 1:v1
            symmetric_gs_u!(u, b)
        end

        # Compute residual
        r = b - actionA_u(u)
        if !quiet
        println("At level $level, After presmoothing,  residual: ", norm(r)/norm(b))  
        println("size(r) = ", size(r))
         

        println("Restricting level $level ...")
        end
        # Restrict residual to coarser grid
        r_coarse = restrict_u_loop(r)

        # Allocate space for the error on the coarser grid
        e_coarse = zeros(size(r_coarse))

        # Recursive call on the coarser grid
        e_solved = vcycle_Au(e_coarse, r_coarse, Int(level/2), v1, v2, quiet)

        # Lift the approximation of the error back to the finer grid
        if !quiet
        println("Lifting level $level ...")
        end
        e_fine = lift_u_loop(e_solved)

        # Correct the approximation on the current grid
        u += e_fine

        resnew = norm(b - actionA_u(u))/norm(b)
        if !quiet
        println("At level $level, After eliminating low frequency error, current residual: ", resnew)
        println("Postsmoothing level $level ...")
        end
        # Post-smoothing with Gauss-Seidel
        for _ in 1:v2
            symmetric_gs_u!(u, b)
        end

        resf = norm(b - actionA_u(u))/norm(b)
        if !quiet
        println("At level $level, After postsmoothing, current residual: ", resf)
        end
        return u
        
    end
end


function vcycle_Av(v, b, level, v1, v2, quiet = true)
    if level == 4
        # Solve exactly at the finest level using ground_truth_4
        if !quiet
            println("Reaching level 4...")
        end
        sol = cg_solver_v(b, zeros(4,3), 300, 1e-10, quiet)
        res = norm(b - actionA_v(sol))/norm(b)
        if !quiet
            println("Level 4, current residual: ", res)

        end
        return sol
    else
        # Pre-smoothing with Gauss-Seidel
        if !quiet
            println("Presmoothing level $level ...")
        end

        for _ in 1:v1
            symmetric_gs_v!(v, b)
        end

        # Compute residual
        r = b - actionA_v(v)
        if !quiet
        println("At level $level, After presmoothing,  residual: ", norm(r)/norm(b))  
         

        println("Restricting level $level ...")
        end
        # Restrict residual to coarser grid
        r_coarse = restrict_v_loop(r)

        # Allocate space for the error on the coarser grid
        e_coarse = zeros(size(r_coarse))

        # Recursive call on the coarser grid
        e_solved = vcycle_Av(e_coarse, r_coarse, Int(level/2), v1, v2, quiet)

        # Lift the approximation of the error back to the finer grid
        if !quiet
        println("Lifting level $level ...")
        end
        e_fine = lift_v_loop(e_solved)

        # Correct the approximation on the current grid
        v += e_fine

        resnew = norm(b - actionA_v(v))/norm(b)
        if !quiet
        println("At level $level, After eliminating low frequency error, current residual: ", resnew)
        println("Postsmoothing level $level ...")
        end
        # Post-smoothing with Gauss-Seidel
        for _ in 1:v2
            symmetric_gs_v!(v, b)
        end

        resf = norm(b - actionA_v(v))/norm(b)
        if !quiet
        println("At level $level, After postsmoothing, current residual: ", resf)
        end
        return v
        
    end
end


function precon_Au(b,tol,v1,v2,maxiter=2,quiet=true)
    n = size(b, 2)
    m = size(b, 1)
    if m != (n-1)
        error("There is a problem with shape of u to be preconditioned")
    end
    u0 = zeros(n-1, n)
    u = zeros(n-1, n)
    for iter in 1:maxiter
        if !quiet
            println("V-cycle iteration ", iter)
        end
        
        u = vcycle_Au(u0, b, n, v1, v2, true)
        u0 = u
        r = b - actionA_u(u)
        res = norm(r)/norm(b)
        if !quiet
            println("V-cycle iteration ", iter, " residual: ", res)
        end
        if res < tol 
            if !quiet
                println("V-cycle converged after ", iter, " iterations")
            end
            return u, iter
        end
        
        
    end
    if !quiet
        println("Warning: V-cycle for Au did not converge")
    end
    return u, maxiter
end

function precon_Av(b,tol,v1,v2,maxiter=2,quiet=true)
    n = size(b, 1)
    m = size(b, 2)
    if m != (n-1)
        error("There is a problem with shape of v to be preconditioned")
    end
    v0 = zeros(n, n-1)
    v = zeros(n, n-1)
    for iter in 1:maxiter
        if !quiet
            println("V-cycle iteration ", iter)
        end
        
        v = vcycle_Av(v0, b, n, v1, v2, true)
        v0 = v
        r = b - actionA_v(v)
        res = norm(r)/norm(b)
        if !quiet
            println("V-cycle iteration ", iter, " residual: ", res)
        end
        if res < tol 
            if !quiet
                println("V-cycle converged after ", iter, " iterations")
            end
            return v, iter
        end
        
        
    end
    if !quiet
        println("Warning: V-cycle for Av did not converge")
    end
    return v, maxiter
end