include("actionB.jl")
include("gs.jl")

function dgs!(u,v,p,f,g,d)
    #d is numerical divergence (negative)
    bp1, bp2 = actionB(p)
    rhs1 = f - bp1
    rhs2 = g - bp2
    gs_update_u!(u,rhs1)
    gs_update_v!(v,rhs2)

    #dimension check
    n1 = size(u,1) + 1
    n2 = size(u,2) 
    n3 = size(v,1)
    n4 = size(v,2) + 1
    n5 = size(p,1)
    n6 = size(p,2)
    n7 = size(d,1)
    n8 = size(d,2)
    # these ni must be equal
    if !(n1==n2==n3==n4==n5==n6==n7==n8)
        error("dgs update: dimension mismatch")
    end

    n = n1
    h = 1.0/n
    r = 0.0
    delta = 0.0
    #update center
    for i in 2:n-1
        for j in 2:n-1
            r = -(u[i,j]-u[i-1,j])/h - (v[i,j]-v[i,j-1])/h - d[i,j]
            delta = r * h / 4
            #update u,v
            u[i-1,j] = u[i-1,j] - delta
            u[i,j] = u[i,j] + delta
            v[i,j-1] = v[i,j-1] - delta
            v[i,j] = v[i,j] + delta
            #update p
            p[i,j] = p[i,j] + r
            p[i+1,j] = p[i+1,j] - r/4
            p[i-1,j] = p[i-1,j] - r/4
            p[i,j+1] = p[i,j+1] - r/4
            p[i,j-1] = p[i,j-1] - r/4
            
        end
    end
    #update upper boundary
    for i in 2:n-1
        j = n
        r = -(u[i,j]-u[i-1,j])/h - (0-v[i,j-1])/h - d[i,j]
        delta = r * h / 3
        #update u,v
        u[i-1,j] = u[i-1,j] - delta
        u[i,j] = u[i,j] + delta
        v[i,j-1] = v[i,j-1] - delta
        #update p
        p[i,j] = p[i,j] + r
        p[i+1,j] = p[i+1,j] - r/3
        p[i-1,j] = p[i-1,j] - r/3
        p[i,j-1] = p[i,j-1] - r/3

    end
    #update lower boundary
    for i in 2:n-1
        j = 1
        r = -(u[i,j]-u[i-1,j])/h - (v[i,j]-0)/h - d[i,j]
        delta = r * h / 3
        #update u,v
        u[i-1,j] = u[i-1,j] - delta
        u[i,j] = u[i,j] + delta
        v[i,j] = v[i,j] + delta
        #update p
        p[i,j] = p[i,j] + r
        p[i+1,j] = p[i+1,j] - r/3
        p[i-1,j] = p[i-1,j] - r/3
        p[i,j+1] = p[i,j+1] - r/3
    end

    #update left boundary
    for j in 2:n-1
        i = 1
        r = -(u[i,j]-0)/h - (v[i,j]-v[i,j-1])/h - d[i,j]
        delta = r * h / 3
        #update u,v
        u[i,j] = u[i,j] + delta
        v[i,j-1] = v[i,j-1] - delta
        v[i,j] = v[i,j] + delta
        #update p
        p[i,j] = p[i,j] + r
        p[i+1,j] = p[i+1,j] - r/3
        p[i,j-1] = p[i,j-1] - r/3
        p[i,j+1] = p[i,j+1] - r/3
    end

    #update right boundary
    for j in 2:n-1
        i = n
        r = -(0-u[i-1,j])/h - (v[i,j]-v[i,j-1])/h - d[i,j]
        delta = r * h / 3
        #update u,v
        u[i-1,j] = u[i-1,j] - delta
        v[i,j-1] = v[i,j-1] - delta
        v[i,j] = v[i,j] + delta
        #update p
        p[i,j] = p[i,j] + r
        p[i-1,j] = p[i-1,j] - r/3
        p[i,j-1] = p[i,j-1] - r/3
        p[i,j+1] = p[i,j+1] - r/3
    end

    #update (1,1)
    i = 1
    j = 1
    r = -(u[i,j]-0)/h - (v[i,j]-0)/h - d[i,j]
    delta = r * h / 2
    #update u,v
    u[i,j] = u[i,j] + delta
    v[i,j] = v[i,j] + delta
    #update p
    p[i,j] = p[i,j] + r
    p[i+1,j] = p[i+1,j] - r/2
    p[i,j+1] = p[i,j+1] - r/2

    #update (1,n)
    i = 1
    j = n
    r = -(u[i,j]-0)/h - (0-v[i,j-1])/h - d[i,j]
    delta = r * h / 2
    #update u,v
    u[i,j] = u[i,j] + delta
    v[i,j-1] = v[i,j-1] - delta
    #update p
    p[i,j] = p[i,j] + r
    p[i+1,j] = p[i+1,j] - r/2
    p[i,j-1] = p[i,j-1] - r/2

    #update (n,1)
    i = n
    j = 1
    r = -(0-u[i-1,j])/h - (v[i,j]-0)/h - d[i,j]
    delta = r * h / 2
    #update u,v
    u[i-1,j] = u[i-1,j] - delta
    v[i,j] = v[i,j] + delta
    #update p
    p[i,j] = p[i,j] + r
    p[i-1,j] = p[i-1,j] - r/2
    p[i,j+1] = p[i,j+1] - r/2

    #update (n,n)
    i = n
    j = n
    r = -(0-u[i-1,j])/h - (0-v[i,j-1])/h - d[i,j]
    delta = r * h / 2
    #update u,v
    u[i-1,j] = u[i-1,j] - delta
    v[i,j-1] = v[i,j-1] - delta
    #update p
    p[i,j] = p[i,j] + r
    p[i-1,j] = p[i-1,j] - r/2
    p[i,j-1] = p[i,j-1] - r/2


end