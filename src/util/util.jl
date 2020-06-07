function Base.convert(T::Type{Vector{SVector{N2,T1}}}, M::Matrix{T2}) where {N2,T1,T2}
    N1 = size(M)[1]
    @assert size(M)[2] == N2
    Mout = [@SVector zeros(T1, N2) for i = 1:N1]
    for i = 1:N1
        Mout[i] = convert(SVector{N2,T1}, M[i,:])
    end
    return Mout
end

function deleteat(M::Array{T,2}, i1::Integer, i2::Integer) where T
    return [M[1:i1 - 1,1:i2 - 1] M[1:i1 - 1,i2 + 1:end];M[i1 + 1:end,1:i2 - 1] M[i1 + 1:end,i2 + 1:end]]
end

deleteat(M::Array{T,2},i::Integer) where T = deleteat(M, i, i)

skew(v::AbstractVector{T}) where T = SMatrix{3,3,T,9}(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0)
skewplusdiag(v::AbstractVector{T},w::T) where T = SMatrix{3,3,T,9}(w, v[3], -v[2], -v[3], w, v[1], v[2], -v[1], w)

Base.:*(u::LinearAlgebra.AdjointAbsVec, v::SVector{1,T}) where T = u * v[1]
Base.:*(u::AbstractVector, v::SVector{1,T}) where T = u * v[1]

@inline svcat(a::T) where T = SVector{1,T}(a)
@inline svcat(a::StaticArray) = a
@inline svcat(a::T, b::T) where T = SVector{2,T}(a, b)
@inline svcat(a::StaticArray, b::StaticArray) = vcat(a, b)
@inline svcat(a::StaticArray{Tuple{N},T,1}, b::T) where {T,N} = vcat(a, SVector{1,T}(b))
@inline svcat(a::T, b::StaticArray{Tuple{N},T,1}) where {T,N} = vcat(SVector{1,T}(a), b)

@inline svcat(a, b, c...) = svcat(svcat(a, b), svcat(c...))

function getfieldnumber(obj)
    i = 1
    while true
        !isdefined(obj, i) ? break : (i+=1)
    end
    return i-1
end

sisnan(a::StaticArray) = any(isnan.(a))

# To fix StaticArray bug
zerodimstaticadjoint(A) = A'
zerodimstaticadjoint(::SMatrix{0,N,T,0}) where {T,N} = SMatrix{N,0,T,0}()
zerodimstaticadjoint(::SMatrix{N,0,T,0}) where {T,N} = SMatrix{0,N,T,0}()