mutable struct Socket{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    cid::Int64

    function Socket(link1::AbstractLink{T},link2::AbstractLink{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 3
        vertices = (p1,p2)
        cid = link2.id

        new{T,Nc}(vertices,cid), link1.id, link2.id
    end
end

@inline function g(joint::Socket,link1::Link,link2::Link,dt,No)
    vertices = joint.vertices
    getx3(link1,dt) + vrotate(vertices[1],getq3(link1,dt)) - (getx3(link2,dt) + vrotate(vertices[2],getq3(link2,dt)))
end

@inline function ∂g∂posa(joint::Socket{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        X = SMatrix{3,3,T,9}(I)

        q = link1.q[No]
        R = 2*VRᵀmat(q)*Rmat(Quaternion(joint.vertices[1]))*LVᵀmat(q)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Socket{T},link1::AbstractLink,link2::Link,No) where T
    if link2.id == joint.cid
        X = SMatrix{3,3,T,9}(-I)

        q = link2.q[No]
        R = -2*VRᵀmat(q)*Rmat(Quaternion(joint.vertices[2]))*LVᵀmat(q)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Socket{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = SMatrix{3,3,T,9}(dt*I)

        q = link1.q[No]
        Ω = 2*dt^2/4*VRᵀmat(q)*Lmat(q)*Rᵀmat(ωbar(link1,dt))*Rmat(Quaternion(joint.vertices[1]))*derivωbar(link1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Socket{T},link1::AbstractLink,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = SMatrix{3,3,T,9}(-dt*I)

        q = link2.q[No]
        Ω = -2*dt^2/4*VRᵀmat(q)*Lmat(q)*Rᵀmat(ωbar(link2,dt))*Rmat(Quaternion(joint.vertices[2]))*derivωbar(link2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function g(joint::Socket,link1::Origin,link2::Link,dt,No)
    vertices = joint.vertices
    vertices[1] - (getx3(link2,dt) + vrotate(vertices[2],getq3(link2,dt)))
end
