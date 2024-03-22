function actionBT(u,v)
    #get the size of the input
    n = size(u, 1) + 1
    m = size(u,2)
    if m != n 
        error("There is a problem with shape of u to be acted by BT")
    end
    n2 = size(v, 1)
    if n2 != n
        error("There is a problem with shape of v to be acted by BT")
    end
    m2 = size(v, 2) + 1
    if m2 != n 
        error("There is a problem with shape of v to be acted by BT")
    end

    h = 1/(n)

    #initialize the output
    w = zeros(n,n)
    # center
    for i in 2:n-1
        for j in 2:n-1
            w[i,j] = (u[i,j]-u[i-1,j]+v[i,j]-v[i,j-1])/h
        end
    end

    # left edge
    for j in 2:n-1
        w[1,j] = (u[1,j]+v[1,j]-v[1,j-1])/h
    end

    # right edge
    for j in 2:n-1
        w[n,j] = (-u[n-1,j]+v[n,j]-v[n,j-1])/h
    end

    # bottom edge
    for i in 2:n-1
        w[i,1] = (u[i,1]-u[i-1,1]+v[i,1])/h
    end

    # top edge
    for i in 2:n-1
        w[i,n] = (u[i,n]- u[i-1,n] - v[i,n-1])/h
    end

    # 4 corners
    w[1,1] = (u[1,1]+v[1,1])/h
    w[1,n] = (u[1,n]-v[1,n-1])/h
    w[n,1] = (-u[n-1,1]+v[n,1])/h
    w[n,n] = (-u[n-1,n]-v[n,n-1])/h

    return -w


end