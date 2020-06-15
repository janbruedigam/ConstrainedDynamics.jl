function Base.convert(T::Type{Vector{SVector{N2,T1}}}, M::Matrix{T2}) where {N2,T1,T2}
    N1 = size(M)[1]
    @assert size(M)[2] == N2
    Mout = [szeros(T1, N2) for i = 1:N1]
    for i = 1:N1
        Mout[i] = convert(SVector{N2,T1}, M[i,:])
    end
    return Mout
end


Base.:*(u::LinearAlgebra.AdjointAbsVec, v::SVector{1,T}) where T = u * v[1]
Base.:*(u::AbstractVector, v::SVector{1,T}) where T = u * v[1]


@inline svcat(a::T) where T = SVector{1,T}(a)
@inline svcat(a::StaticArray) = a
@inline svcat(a::T, b::T) where T = SVector{2,T}(a, b)
@inline svcat(a::StaticArray, b::StaticArray) = vcat(a, b)
@inline svcat(a::StaticArray{Tuple{N},T,1}, b::T) where {T,N} = vcat(a, SVector{1,T}(b))
@inline svcat(a::T, b::StaticArray{Tuple{N},T,1}) where {T,N} = vcat(SVector{1,T}(a), b)
@inline svcat(a, b, c...) = svcat(svcat(a, b), svcat(c...))


@inline szeros(::Type{T}, N) where T = @SVector zeros(T, N)
@inline szeros(N)= @SVector zeros(N)
@inline szeros(::Type{T}, N1, N2) where T = @SMatrix zeros(T, N1, N2)
@inline szeros(N1, N2)= @SMatrix zeros(N1, N2)

@inline sones(::Type{T}, N) where T = @SVector ones(T, N)
@inline sones(N)= @SVector ones(N)
@inline sones(::Type{T}, N1, N2) where T = @SMatrix ones(T, N1, N2)
@inline sones(N1, N2)= @SMatrix ones(N1, N2)

@inline srand(::Type{T}, N) where T = @SVector rand(T, N)
@inline srand(N)= @SVector rand(N)
@inline srand(::Type{T}, N1, N2) where T = @SMatrix rand(T, N1, N2)
@inline srand(N1, N2)= @SMatrix rand(N1, N2)


sisnan(a::StaticArray) = any(isnan.(a))


# To fix StaticArray bug
zerodimstaticadjoint(A) = A'
zerodimstaticadjoint(::SMatrix{0,N,T,0}) where {T,N} = SMatrix{N,0,T,0}()
zerodimstaticadjoint(::SMatrix{N,0,T,0}) where {T,N} = SMatrix{0,N,T,0}()