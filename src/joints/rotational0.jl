mutable struct Rotational0{T,Nc} <: Joint{T,Nc}
    offset::SVector{3,T}
    cid::Int64

    function Rotational0(body1::AbstractBody{T},body2::AbstractBody{T};offset::Quaternion{T}=Quaternion{T}()) where T
        Nc = 3
        offset = Vmat(offset)
        cid = body2.id

        new{T,Nc}(offset,cid), body1.id, body2.id
    end
end

@inline g(joint::Rotational0,body1::Body,body2::Body,dt,No) = (VLᵀmat(getq3(body1,dt))*getq3(body2,dt))-joint.offset

@inline function ∂g∂posa(joint::Rotational0{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T,3,3)

        R = -VRmat(body2.q[No])*RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational0{T},body1::Body,body2::Body,No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T,3,3)

        R = VLᵀmat(body1.q[No])*LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational0{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T,3,3)

        Ω = dt^2/4*VRmat(ωbar(body2,dt))*Rmat(body2.q[No])*Rᵀmat(body1.q[No])*Tmat(T)*derivωbar(body1,dt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational0{T},body1::Body,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T,3,3)

        Ω = dt^2/4*VLᵀmat(ωbar(body1,dt))*Lᵀmat(body1.q[No])*Lmat(body2.q[No])*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline g(joint::Rotational0,body1::Origin,body2::Body,dt,No) = Vmat(getq3(body2,dt))-joint.offset

@inline function ∂g∂posb(joint::Rotational0{T},body1::Origin,body2::Body,No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T,3,3)

        R = VLmat(body2.q[No])*Vᵀmat(T)

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational0{T},body1::Origin,body2::Body,dt,No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T,3,3)

        Ω = dt/2*VLmat(body2.q[No])*derivωbar(body2,dt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
