# No idea what kind of joint this actually is...

mutable struct Rotational2{T,Nc} <: Joint{T,Nc}
    plane::SVector{3,T}
    cid::Int64

    function Rotational2(body1::AbstractBody{T},body2::AbstractBody{T},axis::AbstractVector{T}) where T
        Nc = 1
        cid = body2.id

        new{T,Nc}(axis,cid), body1.id, body2.id
    end
end

@inline g(joint::Rotational2,body1::Body,body2::Body,dt,No) = joint.plane'*(VLᵀmat(getq3(body1,dt))*getq3(body2,dt))

@inline function ∂g∂posa(joint::Rotational2{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T,1,3)

        R = -joint.plane'*VRmat(body2.q[No])*RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational2{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T,1,3)

        R = joint.plane'*VLᵀmat(body1.q[No])*LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational2{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T,1,3)

        Ω = dt^2/4*joint.plane'*VRmat(ωbar(body2,dt))*Rmat(body2.q[No])*Rᵀmat(body1.q[No])*Tmat(T)*derivωbar(body1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational2{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T,1,3)

        Ω = dt^2/4*joint.plane'*VLᵀmat(ωbar(body1,dt))*Lᵀmat(body1.q[No])*Lmat(body2.q[No])*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline g(joint::Rotational2,body1::Origin,body2::Body,dt,No) = joint.plane'*Vmat(getq3(body2,dt))

@inline function ∂g∂posb(joint::Rotational2{T},body1::Origin,body2::Body,No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T,1,3)

        R = joint.plane'*VLmat(body2.q[No])*Vᵀmat(T)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational2{T},body1::Origin,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T,1,3)

        Ω = dt/2*joint.plane'*VLmat(body2.q[No])*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
