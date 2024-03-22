function gs_update_u!(u, b)
    # Get the size of the input
    m = size(u, 1)
    n = size(u, 2)
    
    # Validate the shapes of `u` and `b`
    if m != (n-1) || size(b, 1) != m || size(b, 2) != n
        error("There is a problem with the shape of `u` to be updated by Gauss-Seidel")
    end

    h = 1/n
    # Update the matrix
    for i in 1:n-1
        for j in 1:n
            if i > 1 && i < (n-1) && j > 1 && j < n
                # Center
                u[i,j] = (h^2*b[i,j] + u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1])/4
            elseif i == 1 && j > 1 && j < n
                # Left edge (excluding corners)
                u[i,j] = (h^2*b[i,j] + u[i+1,j] + u[i,j-1] + u[i,j+1])/4
            elseif i == (n-1) && j > 1 && j < n
                # Right edge (excluding corners)
                u[i,j] = (h^2*b[i,j] + u[i-1,j] + u[i,j-1] + u[i,j+1])/4
            elseif j == 1 && i > 1 && i < (n-1)
                # Bottom edge (excluding corners)
                u[i,j] = (h^2*b[i,j] + u[i-1,j] + u[i+1,j] + u[i,j+1])/3
            elseif j == n && i > 1 && i < (n-1)
                # Top edge (excluding corners)
                u[i,j] = (h^2*b[i,j] + u[i-1,j] + u[i+1,j] + u[i,j-1])/3
            elseif (i == 1 && j == 1) 
                # Bottom left corner
                u[i,j] = (h^2*b[i,j] + u[i+1,j] + u[i,j+1])/3
            elseif (i == 1 && j == n)
                # upper left corner
                u[i,j] = (h^2*b[i,j] + u[i+1,j] + u[i,j-1])/3
            elseif (i == (n-1) && j == 1)
                # bottom right corner
                u[i,j] = (h^2*b[i,j] + u[i-1,j] + u[i,j+1])/3
            elseif (i == (n-1) && j == n)
                # upper right corner
                u[i,j] = (h^2*b[i,j] + u[i-1,j] + u[i,j-1])/3
            end

        end
    end
end


function gs_update_v!(v,b)
    # Get the size of the input
    n = size(v, 1)
    m = size(v, 2)
    
    # Validate the shapes of `v` and `b`
    if m != (n-1) || size(b, 1) != n || size(b, 2) != (n-1)
        error("There is a problem with the shape of `v` to be updated by Gauss-Seidel")
    end

    h = 1/n
    # Update the matrix
    for i in 1:n
        for j in 1:n-1
            if i > 1 && i < n && j > 1 && j < (n-1)
                # Center
                v[i,j] = (h^2*b[i,j] + v[i-1,j] + v[i+1,j] + v[i,j-1] + v[i,j+1])/4
            elseif i == 1 && j > 1 && j < (n-1)
                # Left edge (excluding corners)
                v[i,j] = (h^2*b[i,j] + v[i+1,j] + v[i,j-1] + v[i,j+1])/3
            elseif i == n && j > 1 && j < (n-1)
                # Right edge (excluding corners)
                v[i,j] = (h^2*b[i,j] + v[i-1,j] + v[i,j-1] + v[i,j+1])/3
            elseif j == 1 && i > 1 && i < n
                # Bottom edge (excluding corners)
                v[i,j] = (h^2*b[i,j] + v[i-1,j] + v[i+1,j] + v[i,j+1])/4
            elseif j == (n-1) && i > 1 && i < n
                # Top edge (excluding corners)
                v[i,j] = (h^2*b[i,j] + v[i-1,j] + v[i+1,j] + v[i,j-1])/4
            elseif (i == 1 && j == 1) 
                # Bottom left corner
                v[i,j] = (h^2*b[i,j] + v[i+1,j] + v[i,j+1])/3
            elseif (i == 1 && j == (n-1))
                # upper left corner
                v[i,j] = (h^2*b[i,j] + v[i+1,j] + v[i,j-1])/3
            elseif (i == n && j == 1)
                # bottom right corner
                v[i,j] = (h^2*b[i,j] + v[i-1,j] + v[i,j+1])/3
            elseif (i == n && j == (n-1))
                # upper right corner
                v[i,j] = (h^2*b[i,j] + v[i-1,j] + v[i,j-1])/3
            end
        end
    end
end
   