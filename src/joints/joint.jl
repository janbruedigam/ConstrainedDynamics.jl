abstract type Joint{T,Nc} end

Base.show(io::IO, joint::Joint) = summary(io, joint)

getT(joint::Joint{T}) where T = T
getNc(joint::Joint{T,Nc}) where {T,Nc} = Nc

@inline g(joint::Joint{T,Nc}) where {T,Nc} = @SVector zeros(T, Nc)

@inline ∂g∂posa(joint::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T, Nc, 6)
@inline ∂g∂posb(joint::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T, Nc, 6)
@inline ∂g∂vela(joint::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T, Nc, 6)
@inline ∂g∂velb(joint::Joint{T,Nc}) where {T,Nc} = @SMatrix zeros(T, Nc, 6)
