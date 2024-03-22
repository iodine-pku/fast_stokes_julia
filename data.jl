
function gen_rhs(n)
    f = zeros(n-1, n)
    g = zeros(n, n-1)

    h = 1/(n)

    for i in 1:n-1
        for j in 1:n
            f[i,j] = ff(i*h, (j-0.5)*h)
        end
    end

    for i in 1:n-1
        f[i,1] = f[i,1] + fb(i*h)/h
    end

    for i in 1:n-1
        f[i,n] = f[i,n] + ft(i*h)/h
    end

    for i in 1:n
        for j in 1:n-1
            g[i,j] = fg((i-0.5)*h, j*h)
        end
    end

    for j in 1:n-1
        g[1,j] = g[1,j] + fl(j*h)/h
    end

    for j in 1:n-1
        g[n,j] = g[n,j] + fr(j*h)/h
    end

    return f, g

end

function gen_real_sol(n)
    u = zeros(n-1, n)
    v = zeros(n, n-1)
    p = zeros(n, n)
    h = 1/(n)
    for i in 1:n-1
        for j in 1:n
            u[i,j] = fu(i*h, (j-0.5)*h)
        end
    end

    for i in 1:n
        for j in 1:n-1
            v[i,j] = fv((i-0.5)*h, j*h)
        end
    end

    for i in 1:n
        for j in 1:n
            p[i,j] = fp((i-0.5)*h, (j-0.5)*h)
        end
    end

    return u, v, p

end    


function pde_error_sol(u, v)
    # Determine the size of the grid.
    n = size(u, 1) + 1
    h = 1/(n)
    n2 = size(u, 2)
    if n != n2
        error("The size of u has a problem in calculating error.")
    end
    n3 = size(v, 1)
    if n != n3
        error("The size of v has a problem in calculating error..")
    end
    n4 = size(v, 2) + 1
    if n != n4
        error("The size of v has a problem in calculating error..")
    end

    # Generate the real solutions.
    realu, realv, _ = gen_real_sol(n)

    # Calculate the error for u and v.
    error_u = sum((u .- realu).^2)
    error_v = sum((v .- realv).^2)

    # Sum the errors, take the square root, and multiply by h.
    return sqrt(error_u + error_v) * h
end





function ff(x, y)
    return -4π^2 * (2 * cos(2π*x) - 1) * sin(2π*y) + x^2
end

function fg(x, y)
    return 4π^2 * (2 * cos(2π*y) - 1) * sin(2π*x)
end

# Function for expression b
function fb(x)
    return 2π * (cos(2π*x) - 1)
end

# Function for expression t
function ft(x)
    return 4π * sin(π*x)^2
end

# Function for expression l
function fl(y)
    return 4π * sin(π*y)^2
end

# Function for expression r
function fr(y)
    return 2π * (cos(2π*y) - 1)
end

# Function u(x, y)
function fu(x, y)
    return (1 - cos(2π*x)) * sin(2π*y)
end

# Function v(x, y)
function fv(x, y)
    return -(1 - cos(2π*y)) * sin(2π*x)
end

# Function p(x, y)
function fp(x, y)
    return x^3 / 3 - 1 / 12
end
