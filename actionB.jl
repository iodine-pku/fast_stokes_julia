function actionB(p)
    n = size(p, 1)
    m = size(p, 2)
    if n != m
        error("There is a problem with shape of p to be acted by B")
    end

    h = 1/(n)

    #initialize the output
    ans_u = zeros(n-1,n)
    ans_v = zeros(n,n-1)

    for i in 1:n-1
        for j in 1:n
            ans_u[i,j] = (p[i+1,j]-p[i,j])/h
        end
    end

    for i in 1:n
        for j in 1:n-1
            ans_v[i,j] = (p[i,j+1]-p[i,j])/h
        end
    end

    return ans_u, ans_v
end