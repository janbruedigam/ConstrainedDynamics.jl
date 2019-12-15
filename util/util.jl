function Base.convert(T::Type{Vector{SVector{N,T1}}},M::Matrix{T2}) where {N,T1,T2}
    s = size(M)
    @assert s[2] == N
    Mout = [@SVector zeros(T1,N) for i=1:N]
    for i=1:N1
        Mout[i] = convert(SVector{N,T},M[i,:])
    end
    return Mout
end

function Base.deleteat!(M::Array{T,2},i1::Int64,i2::Int64) where T
    M[1:end-1,1:end-1] = [M[1:i1-1,1:i2-1] M[1:i1-1,i2+1:end];M[i1+1:end,1:i2-1] M[i1+1:end,i2+1:end]]
end

Base.deleteat!(M::Array,i::Int64) = deleteat!(M,i,i)
