mutable struct Translational2{T,Nc} <: Joint{T,Nc}
    vertices::NTuple{2,SVector{3,T}}
    plane::SVector{3,T}
    cid::Int64

    function Translational2(body1::AbstractBody{T},body2::AbstractBody{T},p1::AbstractVector{T},p2::AbstractVector{T},axis::AbstractVector{T}) where T
        Nc = 1
        vertices = (p1,p2)
        cid = body2.id

        new{T,Nc}(vertices,axis,cid), body1.id, body2.id
    end
end

@inline function g(joint::Translational2,body1::Body,body2::Body,dt,No)
    vertices = joint.vertices
    q1 = getq3(body1,dt)
    joint.plane'*vrotate(getx3(body2,dt) + vrotate(vertices[2],getq3(body2,dt)) - (getx3(body1,dt) + vrotate(vertices[1],q1)),inv(q1))
end

@inline function ∂g∂posa(joint::Translational2{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        point2 = body2.x[No] + vrotate(joint.vertices[2],body2.q[No])

        X = -joint.plane'*VLᵀmat(q1)*RVᵀmat(q1)

        R = joint.plane'*2*VLᵀmat(q1)*(Lmat(Quaternion(point2)) - Lmat(Quaternion(body1.x[No])))*LVᵀmat(q1)

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational2{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        q2 = body2.q[No]

        X = joint.plane'*VLᵀmat(q1)RVᵀmat(q1)

        R = joint.plane'*2*VLᵀmat(q1)*Rmat(q1)*Rᵀmat(q2)*Rmat(Quaternion(joint.vertices[2]))*LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational2{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        ωbar1 = ωbar(body1,dt)
        point2 = body2.x[No] + dt*getvnew(body2) + dt^2/4*vrotate(vrotate(joint.vertices[2],ωbar(body2,dt)),body2.q[No])

        V = -dt^3/4*joint.plane'*VLᵀmat(ωbar1)Lᵀmat(q1)Rmat(ωbar1)RVᵀmat(q1)

        Ω = 2*dt^2/4*joint.plane'*VLᵀmat(ωbar1)*Lᵀmat(q1)*(Lmat(Quaternion(point2))-Lmat(Quaternion(body1.x[No] + dt*getvnew(body1))))*Lmat(q1)*derivωbar(body1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational2{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        q1 = body1.q[No]
        q2 = body2.q[No]
        ωbar1 = ωbar(body1,dt)
        V = dt^3/4*joint.plane'*VLᵀmat(ωbar1)Lᵀmat(q1)Rmat(ωbar1)RVᵀmat(q1)

        Ω = 2*dt^4/16*joint.plane'*VLᵀmat(ωbar1)*Lᵀmat(q1)*Lmat(q2)*Rmat(ωbar1)*Rmat(q1)*Rᵀmat(q2)*Rᵀmat(ωbar(body2,dt))*Rmat(Quaternion(joint.vertices[2]))*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function g(joint::Translational2,body1::Origin,body2::Body,dt,No)
    vertices = joint.vertices
    joint.plane'*(getx3(body2,dt) + vrotate(vertices[2],getq3(body2,dt)) - vertices[1])
end

@inline function ∂g∂posb(joint::Translational2{T},body1::Origin,body2::Body,No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        X = joint.plane'

        R = joint.plane'*2*VRᵀmat(q2)*Rmat(Quaternion(joint.vertices[2]))*LVᵀmat(q2)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Translational2{T},body1::Origin,body2::Body,dt,No) where T
    if body2.id == joint.cid
        q2 = body2.q[No]

        V = dt*joint.plane'

        Ω = 2*dt^2/4*joint.plane'*VLmat(q2)*Rᵀmat(q2)*Rᵀmat(ωbar(body2,dt))*Rmat(Quaternion(joint.vertices[2]))*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
