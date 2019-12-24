mutable struct SocketYZ{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    cid::Int64

    function SocketYZ(link1::AbstractLink{T},link2::AbstractLink{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 2
        vertices = (p1,p2)
        cid = link2.id

        new{T,Nc}(vertices,cid), link1.id, link2.id
    end
end

@inline function g(joint::SocketYZ,link1::Link,link2::Link,dt,No)
    vertices = joint.vertices
    (getx3(link1,dt) + vrotate(vertices[1],getq3(link1,dt)) - (getx3(link2,dt) + vrotate(vertices[2],getq3(link2,dt))))[SVector{2,Int64}(2,3)]
end

@inline function ∂g∂posa(joint::SocketYZ{T},link1::Link,link2::Link,No) where T
    if link2.id == joint.cid
        X = SMatrix{2,3,T,6}(0,0, 1,0, 0,1)

        q = link1.q[No]
        R = 2*(VRᵀmat(q)*Rmat(Quaternion(joint.vertices[1]))*LVᵀmat(q))[SVector{2,Int64}(2,3),:]

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::SocketYZ{T},link1::AbstractLink,link2::Link,No) where T
    if link2.id == joint.cid
        X = SMatrix{2,3,T,6}(0,0, -1,0, 0,-1)

        q = link2.q[No]
        R = -2*(VRᵀmat(q)*Rmat(Quaternion(joint.vertices[2]))*LVᵀmat(q))[SVector{2,Int64}(2,3),:]

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::SocketYZ{T},link1::Link,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = SMatrix{2,3,T,6}(0,0, dt,0, 0,dt)

        q = link1.q[No]
        Ω = 2*dt^2/4*(VRᵀmat(q)*Lmat(q)*Rᵀmat(ωbar(link1,dt))*Rmat(Quaternion(joint.vertices[1]))*derivωbar(link1,dt))[SVector{2,Int64}(2,3),:]

        return [V Ω]
    else
        return ∂g∂vel(joint)
    end
end

@inline function ∂g∂velb(joint::SocketYZ{T},link1::AbstractLink,link2::Link,dt,No) where T
    if link2.id == joint.cid
        V = SMatrix{2,3,T,6}(0,0, -dt,0, 0,-dt)

        q = link2.q[No]
        Ω = -2*dt^2/4*(VRᵀmat(q)*Lmat(q)*Rᵀmat(ωbar(link2,dt))*Rmat(Quaternion(joint.vertices[2]))*derivωbar(link2,dt))[SVector{2,Int64}(2,3),:]

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function g(joint::SocketYZ,link1::Origin,link2::Link,dt,No)
    vertices = joint.vertices
    (vertices[1] - (getx3(link2,dt) + vrotate(vertices[2],getq3(link2,dt))))[SVector{2,Int64}(2,3)]
end
