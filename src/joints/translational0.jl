mutable struct Translational0{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    cid::Int64

    function Translational0(body1::AbstractBody{T},body2::AbstractBody{T},p1::AbstractVector{T},p2::AbstractVector{T}) where T
        Nc = 3
        vertices = (p1,p2)
        cid = body2.id

        new{T,Nc}(vertices,cid), body1.id, body2.id
    end
end

@inline function g(joint::Translational0,body1::Body,body2::Body,dt,No)
    vertices = joint.vertices
    getx3(body2,dt) + vrotate(vertices[2],getq3(body2,dt)) - (getx3(body1,dt) + vrotate(vertices[1],getq3(body1,dt)))
end

@inline function ∂g∂posa(joint::Translational0{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        X = SMatrix{3,3,T,9}(-I)

        q = body1.q[No]
        R = -2*VRᵀmat(q)*Rmat(Quaternion(joint.vertices[1]))*LVᵀmat(q)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational0{T},body1::AbstractBody,body2::Body,No) where T
    if body2.id == joint.cid
        X = SMatrix{3,3,T,9}(I)

        q = body2.q[No]
        R = 2*VRᵀmat(q)*Rmat(Quaternion(joint.vertices[2]))*LVᵀmat(q)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational0{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = SMatrix{3,3,T,9}(-dt*I)

        q = body1.q[No]
        Ω = -2*dt^2/4*VRᵀmat(q)*Lmat(q)*Rᵀmat(ωbar(body1,dt))*Rmat(Quaternion(joint.vertices[1]))*derivωbar(body1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational0{T},body1::AbstractBody,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = SMatrix{3,3,T,9}(dt*I)

        q = body2.q[No]
        Ω = 2*dt^2/4*VRᵀmat(q)*Lmat(q)*Rᵀmat(ωbar(body2,dt))*Rmat(Quaternion(joint.vertices[2]))*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function g(joint::Translational0,body1::Origin,body2::Body,dt,No)
    vertices = joint.vertices
    getx3(body2,dt) + vrotate(vertices[2],getq3(body2,dt)) - vertices[1]
end
