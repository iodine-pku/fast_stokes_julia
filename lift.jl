function lift_p_loop(q)
    n = size(q, 1)
    m = size(q,2)

    if m != n
        error("p to be lifted must be square")
    end

    p = zeros(2n, 2n)
    for i in 1:n
        for j in 1:n
            p[2i-1, 2j-1] = q[i, j]
            p[2i-1, 2j] = q[i, j]
            p[2i, 2j-1] = q[i, j]
            p[2i, 2j] = q[i, j]
        end
    end
    return p
end

function lift_p_ind(q)
    n = size(q, 1)
    m = size(q,2)

    if m != n
        error("Matrix must be square")
    end

    p = zeros(2n, 2n)
    p[1:2:end-1, 1:2:end-1] .= q
    p[1:2:end-1, 2:2:end] .= q
    p[2:2:end, 1:2:end-1] .= q
    p[2:2:end, 2:2:end] .= q
    return p
end

function lift_u_loop(q)
    n = size(q, 2)
    m = size(q,1)

    if m != (n-1)
        error("There is a problem with shape of u to be lifted")
    end

    u = zeros(2n-1, 2n)
    for j in 1:n
        for i in 1:n-1
            u[2i,2j-1] = q[i,j]
            u[2i,2j] = q[i,j]
        end    
    end

    for j in 1:n
        u[1,2j-1] = 0.5 * q[1,j]
        u[1,2j] = 0.5 * q[1,j]
        u[2n-1,2j-1] = 0.5 * q[n-1,j]
        u[2n-1,2j] = 0.5 * q[n-1,j]
    end

    for j in 1:n
        for i in 2:n-1
            u[2i-1,2j-1] = 0.5 * (q[i,j]+q[i-1,j])
            u[2i-1,2j] = 0.5 * (q[i,j]+q[i-1,j])
            
        end
    end
    return u
end

function lift_v_loop(q)
    n = size(q, 1)
    m = size(q,2)

    if m != (n-1)
        error("There is a problem with shape of v to be lifted")
    end

    v = zeros(2n, 2n-1)
    for i in 1:n
        for j in 1:n-1
            v[2i-1,2j] = q[i,j]
            v[2i,2j] = q[i,j]
        end    
    end

    for i in 1:n
        v[2i-1,1] = 0.5 * q[i,1]
        v[2i,1] = 0.5 * q[i,1]
        v[2i-1,2n-1] = 0.5 * q[i,n-1]
        v[2i,2n-1] = 0.5 * q[i,n-1]
    end

    for i in 1:n
        for j in 2:n-1
            v[2i-1,2j-1] = 0.5 * (q[i,j]+q[i,j-1])
            v[2i,2j-1] = 0.5 * (q[i,j]+q[i,j-1])
        end
    end
    return v
end
