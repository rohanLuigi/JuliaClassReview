function Jacobi(A,x,b,p;tol=0.000001,output=10)
    D       = diag(A)
    res     = A*x
    resid   = b-res
    tempSol = resid./D
    diff    = sum(resid.^2)
    n    = size(A,1)
    iter = 0
    while ((diff/n > tol) & (iter<1000))
        iter += 1
        x = p*tempSol + (1-p)*x
        res     = A*x
        resid   = b-res
        tempSol = resid./D + x
        diff    = sum(resid.^2)
        if iter%output == 0
            println(iter," ",diff/n)
        end
    end
    return x
end

function GaussSeidel(A,x,b;tol=0.000001,output=10)
    n = size(x,1)
    iter = 0
    diff = 1.0
    while ((diff/n > tol) & (iter<1000))
        iter += 1
        for i=1:n
            x[i] = ((b[i] - A[i,:]*x)/A[i,i])[1,1] + x[i]
        end
        diff = sum((A*x-b).^2)
        println(iter," ",diff/n)
    end
    return x
end

function GaussSeidel(A,x,b;tol=0.000001,output=10)
    n = size(x,1)
    for i=1:n
        x[i] = ((b[i] - A[i,:]*x)/A[i,i])[1,1] + x[i]
    end
    diff = sum((A*x-b).^2)
    iter = 0
    while ((diff/n > tol) & (iter<1000))
        iter += 1
        for i=1:n
            x[i] = ((b[i] - A[i,:]*x)/A[i,i])[1,1] + x[i]
        end
        diff = sum((A*x-b).^2)
        if iter%output == 0
            println(iter," ",diff/n)
        end
    end
    return x
end

