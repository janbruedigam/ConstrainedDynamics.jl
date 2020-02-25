abstract type Contact{T} end

# Base.show(io::IO, joint::Joint) = summary(io, joint)

getT(contact::Contact{T}) where T = T

# @inline g(joint::Joint{T,Nc}) where {T,Nc} = @SVector zeros(T,Nc)
#
# @inline ∂g∂posa(contact::Contact{T}) where {T} = @SMatrix zeros(T,Nc,6)
# @inline ∂g∂posb(contact::Contact{T}) where {T} = @SMatrix zeros(T,Nc,6)
# @inline ∂g∂vela(contact::Contact{T}) where {T} = @SMatrix zeros(T,Nc,6)
# @inline ∂g∂velb(contact::Contact{T}) where {T} = @SMatrix zeros(T,Nc,6)
