function restrict_p_loop(p)
    n = size(p, 1) 
    m = size(p, 2) 


    if m != n
        error("p to be restricted must be square")
    end
    if n % 2 != 0
        error("p to be restricted must have an even number of rows and columns")
    end
    n = n ÷ 2

    q = zeros(n, n)
    for i in 1:n
        for j in 1:n
            q[i, j] = 0.25 * (p[2*i-1, 2*j-1] + p[2*i-1, 2*j] + p[2*i, 2*j-1] + p[2*i, 2*j])
        end
    end
    return q
end

function restrict_p_ind(p)
    n = size(p, 1) ÷ 2
    m = size(p, 2) ÷ 2
    if m != n
        error("p to be restricted must be square")
    end
    q = 0.25 .* (
        p[1:2:end-1, 1:2:end-1] .+ p[1:2:end-1, 2:2:end] .+
        p[2:2:end, 1:2:end-1] .+ p[2:2:end, 2:2:end]
    )
    return q
end

function restrict_u_loop(u)
    n = size(u,2)
    m = size(u,1)
    if m != (n-1)
        error("There is a problem with shape of u to be restricted")
    end
    if n % 2 != 0
        error("u to be restricted must have an even number of columns")
    end
    n = n ÷ 2


    q = zeros(n-1, n)
    for i in 1:(n-1)
        for j in 1:n
            q[i, j] = 0.25 * (u[2*i, 2*j-1] + u[2*i, 2*j]) +
                      0.125 * (u[2*i-1, 2*j-1] + u[2*i+1, 2*j-1] +
                               u[2*i-1, 2*j] + u[2*i+1, 2*j])
        end
    end
    return q
end

function restrict_v_loop(v)
    n = size(v, 1) 
    m = size(v, 2) 
    if m != (n-1)
        error("There is a problem with shape of v to be restricted")
    end

    if n % 2 != 0
        error("v to be restricted must have an even number of rows")
    end
    n = n ÷ 2
    
    q = zeros(n, n-1)
    for i in 1:n
        for j in 1:(n-1)
            q[i, j] = 0.25 * (v[2*i-1, 2*j] + v[2*i, 2*j]) +
                      0.125 * (v[2*i-1, 2*j-1] + v[2*i-1, 2*j+1] +
                               v[2*i, 2*j-1] + v[2*i, 2*j+1])
        end
    end
    return q
end
