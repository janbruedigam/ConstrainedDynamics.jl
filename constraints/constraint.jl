using StaticArrays

include("../util/quaternion.jl")
include("../node.jl")
include("../link.jl")

abstract type Constraint{T,Nc,N,Nc²,NcN,Nl} <: Node{T,Nc,N,Nc²,NcN} end

include("fixedposition.jl")
include("fixedorientation.jl")
include("socket.jl")
include("axis.jl")
include("combined.jl")


function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, C::Constraint)
    summary(io, C); println(io)
    print(io, "\nConnected links: ")
    show(io, mime, linkids(C))
end

Base.show(io::IO, C::Link) = summary(io, C)

function Base.getproperty(C::Constraint{T,Nc},x::Symbol) where {T,Nc}
    if x == :T
        return T
    elseif x == :Nc
        return Nc
    else
        return getfield(C,x)
    end
end

@inline g(C::Constraint{T,Nc}) where {T,Nc} = @SVector zeros(T,Nc)

@inline zeroBlock(C::Constraint{T,Nc}) where {T,Nc} = @SMatrix zeros(T,Nc,6)
@inline ∂g∂posa(C::Constraint) = zeroBlock(C)
@inline ∂g∂posb(C::Constraint) = zeroBlock(C)
@inline ∂g∂vela(C::Constraint) = zeroBlock(C)
@inline ∂g∂velb(C::Constraint) = zeroBlock(C)

@inline linkids(C::Constraint{T,Nc,N,Nc²,NcN,1}) where {T,Nc,N,Nc²,NcN} = @SVector [C.link.data.id]
@inline linkids(C::Constraint{T,Nc,N,Nc²,NcN,2}) where {T,Nc,N,Nc²,NcN} = @SVector [C.link1.data.id, C.link2.data.id]

@inline Gtλ(L::Link,C::Constraint) = L.data.id==linkids(C)[1] ? ∂g∂posa(C)'*C.data.s1 : ∂g∂posb(C)'*C.data.s1
