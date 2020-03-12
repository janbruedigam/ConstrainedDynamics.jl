mutable struct Rotational3{T,Nc} <: Joint{T,Nc}
    qoff::Quaternion{T}
    cid::Int64

    function Rotational3(body1::AbstractBody{T}, body2::AbstractBody{T};offset::Quaternion{T} = Quaternion{T}()) where T
        Nc = 3
        cid = body2.id

        new{T,Nc}(offset, cid), body1.id, body2.id
    end
end

function setForce!(joint::Rotational3, body1::AbstractBody, body2::Body, τ, No)
    return
end

@inline minimalCoordinates(joint::Rotational3, body1::Body{T}, body2::Body, No) where T = SVector{0,T}()

# @inline g(joint::Rotational3, body1::Body, body2::Body, Δt, No) = VLᵀmat(getq3(body1, Δt)) * getq3(body2, Δt) - joint.offset
@inline g(joint::Rotational3, body1::Body, body2::Body, Δt, No) = VRmat(joint.qoff) * Lᵀmat(getq3(body1, Δt)) * getq3(body2, Δt)

@inline function ∂g∂posa(joint::Rotational3{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 3, 3)

        # R = -VRmat(body2.q[No]) * RᵀVᵀmat(body1.q[No])
        R = -VRmat(joint.qoff) * Rmat(body2.q[No]) * RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational3{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 3, 3)

        # R = VLᵀmat(body1.q[No]) * LVᵀmat(body2.q[No])
        R = VRmat(joint.qoff) * Lᵀmat(body1.q[No]) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 3, 3)

        # Ω = VRmat(ωbar(body2, Δt)) * Rmat(body2.q[No]) * Rᵀmat(body1.q[No]) * Tmat(T) * derivωbar(body1, Δt)
        Ω = VRmat(joint.qoff) * Rmat(ωbar(body2, Δt)) * Rmat(body2.q[No]) * Rᵀmat(body1.q[No]) * Tmat(T) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 3, 3)

        # Ω = VLᵀmat(ωbar(body1, Δt)) * Lᵀmat(body1.q[No]) * Lmat(body2.q[No]) * derivωbar(body2, Δt)
        Ω = VRmat(joint.qoff) * Lᵀmat(ωbar(body1, Δt)) * Lᵀmat(body1.q[No]) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end

@inline minimalCoordinates(joint::Rotational3, body1::Origin{T}, body2::Body, No) where T = SVector{0,T}()

# @inline g(joint::Rotational3, body1::Origin, body2::Body, Δt, No) = Vmat(getq3(body2, Δt)) - joint.offset
@inline g(joint::Rotational3, body1::Origin, body2::Body, Δt, No) = VRmat(joint.qoff) * getq3(body2, Δt)

@inline function ∂g∂posb(joint::Rotational3{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 3, 3)

        # R = VLmat(body2.q[No]) * Vᵀmat(T)
        R = VRmat(joint.qoff) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational3{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 3, 3)

        # Ω = VLmat(body2.q[No]) * derivωbar(body2, Δt)
        Ω = VRmat(joint.qoff) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
