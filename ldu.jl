using LinearAlgebra

struct LDU{T} <: Factorization{T}
    n::Int64
    A::Array{T,2}
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::LDU)
    summary(io, F); println(io)
    println(io, "L factor:")
    show(io, mime, F.L)
    println(io, "\nD factor:")
    show(io, mime, F.D)
    println(io, "\nU factor:")
    show(io, mime, F.U)
end

function ldu!(A)
    s = size(A)[1]

    for n = 1:s
        # L and R
        for p = 1:n-1
            for k = 1:p-1
                A[n,p] -= A[n,k]*A[k,k]*A[k,p]
                A[p,n] -= A[p,k]*A[k,k]*A[k,n]
            end
            A[n,p] = A[n,p]/A[p,p]
            A[p,n] = A[p,p]\A[p,n]
        end

        # D
        for k = 1:n-1
            A[n,n] -= A[n,k]*A[k,k]*A[k,n]
        end
    end
    LDU(s,A)
end

ldu(A) = ldu!(copy(A))

function Base.getproperty(F::LDU{T},x::Symbol) where T
    if x == :L
        A = F.A
        n = F.n
        L = Matrix{T}(I,n,n)
        for i = 1:n
            for j = 1:i-1
                L[i,j] = A[i,j]
            end
        end
        return L
    elseif x == :D
        D = diagm(diag(F.A))
        return D
    elseif x == :U
        A = F.A
        n = F.n
        R = Matrix{T}(I,n,n)
        for j = 1:n
            for i = 1:j-1
                R[i,j] = A[i,j]
            end
        end
        return R
    else
        return getfield(F,x)
    end
end

function Base.:\(F::LDU,b::AbstractVector)
    s = F.n
    A = F.A
    x = copy(b)

    # L
    for n = 1:s
        for k = 1:n-1
            x[n] -= A[n,k]*x[k]
        end
    end

    # D and R
    for n = s:-1:1
        x[n] /= A[n,n]
        for k = n+1:s
            x[n] -= A[n,k]*x[k]
        end
    end

    return x
end
