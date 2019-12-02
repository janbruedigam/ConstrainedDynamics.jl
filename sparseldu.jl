using LinearAlgebra
using StaticArrays

include("link.jl")
include("constraints/constraint.jl")

@inline function setD!(link::Link{T}) where T
    dynT,dynR = ∂dyn∂vel(link)
    Z = @SMatrix zeros(T,3,3)

    link.data.D = [[dynT Z];[Z dynR]]
    return nothing
end

@inline setD!(C::Constraint{T,Nc}) where {T,Nc} = (C.data.D = @SMatrix zeros(T,Nc,Nc); nothing)
@inline updateD!(node,child) = (d = child.data; node.data.D -= d.JL*d.D*d.JU; nothing)
@inline invertD!(node) = (d = node.data; d.Dinv = inv(d.D); nothing)

@inline function setJ!(L::Link,C::Constraint)
    data = L.data

    if data.id==linkids(C)[1]
        data.JL = ∂g∂vela(C)
        data.JU = -∂g∂posa(C)'
    else
        data.JL = ∂g∂velb(C)
        data.JU = -∂g∂posb(C)'
    end

    return nothing
end

@inline function setJ!(C::Constraint,L::Link)
    data = C.data

    if L.data.id==linkids(C)[1]
        data.JL = -∂g∂posa(C)'
        data.JU = ∂g∂vela(C)
    else
        data.JL = -∂g∂posb(C)'
        data.JU = ∂g∂velb(C)
    end

    return nothing
end

@inline function updateJ!(node)
    d = node.data
    d.JL = d.JL*d.Dinv
    d.JU = d.Dinv*d.JU
    return nothing
end

@inline setSol!(link::Link) = (link.data.ŝ = dynamics(link); nothing)
@inline setSol!(C::Constraint) = (C.data.ŝ = g(C); nothing)

# (A) For extended equations
# @inline addGtλ!(L::Link,C::Constraint) = (L.data.ŝ -= Gtλ(L,C); nothing)
@inline addλ0!(C::Constraint) = (C.data.ŝ += C.data.s0; nothing)

@inline LSol!(node,child) = (d = child.data; node.data.ŝ -= d.JL*d.ŝ; nothing)
@inline DSol!(node) = (d = node.data; d.ŝ = d.Dinv*d.ŝ; nothing)
@inline USol!(node,parent) = (d = node.data; d.ŝ -= d.JU*parent.data.ŝ; nothing)




@inline function setNormf!(link::Link,robot::Robot)
    data = link.data
    data.f = dynamics(link)
    for (i,connected) in enumerate(robot.adjacency[data.id])
        connected ? GtλTof!(link,robot.nodes[i]) : nothing
    end
    data.normf = data.f'*data.f
    return nothing
end

@inline function setNormf!(C::Constraint,::Robot)
    data = C.data
    data.f = g(C)
    data.normf = data.f'*data.f
    return nothing
end

@inline GtλTof!(L::Link,C::Constraint) = (L.data.f -= Gtλ(L,C); nothing)
