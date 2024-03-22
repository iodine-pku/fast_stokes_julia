function cal_res(u,v,p,f,g,d)
    bp1,bp2 = actionB(p)
    res1 = f - actionA_u(u) - bp1
    res2 = g - actionA_v(v) - bp2
    res3 = d - actionBT(u,v)
    #calculate norms
    norm1 =sqrt(norm(res1)^2 + norm(res2)^2+norm(res3)^2)
    norm2 = sqrt(norm(f)^2 + norm(g)^2 + norm(d)^2)
    rela_norm = norm1/norm2


    return res1,res2,res3, rela_norm
end


function general_uzawa_update(p,f,g,d,alpha=1.0)
    n = size(p, 1)
    pnew = zeros(n,n)
    bp1 , bp2 = actionB(p)
    r1 = zeros(n-1,n)
    r2 = zeros(n,n-1)
    r1 = f - bp1
    r2 = g - bp2
    unew = cg_solver_u(r1,zeros(n-1,n))
    vnew = cg_solver_v(r2,zeros(n,n-1))  
  

    pnew = p + alpha*(actionBT(unew,vnew)-d)

    return unew, vnew, pnew

end



function basesolver(f,g,d, tol=1e-8, max_iter=20,silent=true)
    m = size(d,1)
    u = zeros(m-1,m)
    v = zeros(m,m-1)
    p = zeros(m,m)
    unew = zeros(m-1,m)
    vnew = zeros(m,m-1)
    pnew = zeros(m,m)
    
    for i = 1:max_iter
        unew, vnew, pnew = general_uzawa_update(p,f,g,d)
        res1, res2, res3, relae = cal_res(unew,vnew,pnew,f,g,d)
        if !silent
        println("Iteration: ", i, " Residual: ", relae)
        end
        if relae < tol
            if !silent
            println("Converged in ", i, " iterations")
            end
            
            return unew, vnew, pnew, i
        end
        u = unew
        v = vnew
        p = pnew
    end
    if !silent
    println("Basesolver Did not converge")
    end
    return  u, v, p, max_iter
end