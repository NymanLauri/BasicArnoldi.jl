module BasicArnoldi
export basic_arnoldi
    """
    Computes the basic Arnoldi algorithm for a matrix A and a vector v1.
    The parameter m specifies the maximum order of the computed Krylov subspace.
    Stores the Krylov subspace basis vectors in v and the Hessenberg matrix in H.
    """
    function basic_arnoldi(A,v1,m)
        n = length(v1)
        v = Matrix{Float64}(n,m+1)
        H = zeros(m+1,m)
        v[:,1] .= v1 ./ norm(v1)
        for j=1:m
            v_next = view(v, :, j + 1)
            A_mul_B!(v_next, A, view(v, :, j))
            for i=1:j
                H[i,j] = dot(v_next,view(v, :, i))
                v_next .-= H[i,j] .* view(v, :,i)
            end
            H[j+1,j]=norm(v_next)
            v[:,j+1] .= v_next ./ H[j+1,j]
        end
        v,H
    end
end