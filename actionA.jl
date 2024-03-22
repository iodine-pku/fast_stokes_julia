function actionA_u(u)
    #get the size of the input
    n = size(u, 1) + 1
    m = size(u, 2)
    if m != n 
        error("There is a problem with shape of u to be acted by A")
    end
    
    h = 1/(n)
    #initialize the output
    w = zeros(n-1,n)
    # center
    for i in 2:n-2
        for j in 2:n-1
            w[i,j] = (4*u[i,j] - u[i-1,j] - u[i+1,j] - u[i,j-1] - u[i,j+1])/(h^2)
        end
    end
    # left edge
    for j in 2:n-1
        w[1,j] = (4*u[1,j] - u[2,j] - u[1,j-1] - u[1,j+1])/(h^2)
    end

    # right edge
    for j in 2:n-1
        w[n-1,j] = (4*u[n-1,j] - u[n-2,j] - u[n-1,j-1] - u[n-1,j+1])/(h^2)
    end

    # bottom edge
    for i in 2:n-2
        w[i,1] = (3*u[i,1] - u[i-1,1] - u[i+1,1] - u[i,2])/(h^2)
    end

    # top edge
    for i in 2:n-2
        w[i,n] = (3*u[i,n] - u[i-1,n] - u[i+1,n] - u[i,n-1])/(h^2)
    end

    # 4 corners
    w[1,1] = (3*u[1,1] - u[2,1] - u[1,2])/(h^2)
    w[1,n] = (3*u[1,n] - u[2,n] - u[1,n-1])/(h^2)
    w[n-1,1] = (3*u[n-1,1] - u[n-2,1] - u[n-1,2])/(h^2)
    w[n-1,n] = (3*u[n-1,n] - u[n-2,n] - u[n-1,n-1])/(h^2)
    return w

end

function actionA_v(v)
    #get the size of the input
    n = size(v, 1)
    m = size(v, 2) + 1
    if n != m 
        error("There is a problem with shape of v to be acted by A")
    end
    
    h = 1/(n)
    #initialize the output
    w = zeros(n,n-1)
    # center
    for i in 2:n-1
        for j in 2:n-2
            w[i,j] = (4*v[i,j] - v[i-1,j] - v[i+1,j] - v[i,j-1] - v[i,j+1])/(h^2)
        end
    end
    # left edge
    for j in 2:n-2
        w[1,j] = (3*v[1,j] - v[1,j-1] - v[1,j+1] - v[2,j])/(h^2)
    end

    # right edge
    for j in 2:n-2
        w[n,j] = (3*v[n,j] - v[n,j-1] - v[n,j+1] - v[n-1,j])/(h^2)
    end

    # bottom edge
    for i in 2:n-1
        w[i,1] = (4*v[i,1] - v[i-1,1] - v[i+1,1] - v[i,2])/(h^2)
    end

    # top edge
    for i in 2:n-1
        w[i,n-1] = (4*v[i,n-1] - v[i-1,n-1] - v[i+1,n-1] - v[i,n-2])/(h^2)
    end

    # 4 corners
    w[1,1] = (3*v[1,1] - v[1,2] - v[2,1])/(h^2)
    w[1,n-1] = (3*v[1,n-1] - v[1,n-2] - v[2,n-1])/(h^2)
    w[n,1] = (3*v[n,1] - v[n,2] - v[n-1,1])/(h^2)
    w[n,n-1] = (3*v[n,n-1] - v[n,n-2] - v[n-1,n-1])/(h^2)
    return w
end