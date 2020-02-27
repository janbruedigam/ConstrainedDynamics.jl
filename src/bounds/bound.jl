abstract type Bound{T} end

Base.show(io::IO, bound::Bound) = summary(io, bound)

getT(bound::Bound{T}) where T = T


@inline g(bound::Bound{T}) where T = zero(T)

∂g∂pos(bound::Bound{T}) where T = @SVector zeros(T,6)
∂g∂vel(bound::Bound{T}) where T = @SVector zeros(T,6)
schurf(bound::Bound{T}) where T = @SVector zeros(T,6)
schurD(bound::Bound{T}) where T = @SMatrix zeros(T,6,6)

