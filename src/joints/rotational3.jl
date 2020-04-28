@inline function getVelocityDelta(joint::Rotational3, body1::AbstractBody, body2::Body{T}, ω::Union{Nothing,SVector{0,T}}) where T
    Δω = @SVector zeros(T,3)
    return Δω
end

@inline function getPositionDelta(joint::Rotational3, body1::AbstractBody, body2::Body{T}, θ::Union{Nothing,SVector{0,T}}) where T
    Δq = Quaternion{T}()
    return Δq
end

@inline function minimalCoordinates(joint::Rotational3, body1::AbstractBody{T}, body2::Body, No) where T
    SVector{0,T}()
end


@inline function g(joint::Rotational3, body1::Body, body2::Body, Δt, No)
    VLᵀmat(joint.qoff) * Lᵀmat(getq3(body1, Δt)) * getq3(body2, Δt)
end

@inline function g(joint::Rotational3, body1::Origin, body2::Body, Δt, No)
    VLᵀmat(joint.qoff) * getq3(body2, Δt)
end


@inline function ∂g∂posa(joint::Rotational3{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 3, 3)
        R = -VLᵀmat(joint.qoff) * Rmat(body2.q[No]) * RᵀVᵀmat(body1.q[No])

        return [X R]
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Rotational3{T}, body1::Body, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 3, 3)
        R = VLᵀmat(joint.qoff) * Lᵀmat(body1.q[No]) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 3, 3)
        Ω = VLᵀmat(joint.qoff) * Rmat(ωbar(body2, Δt)) * Rmat(body2.q[No]) * Rᵀmat(body1.q[No]) * Tmat(T) * derivωbar(body1, Δt)

        return [V Ω]
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational3{T}, body1::Body, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 3, 3)
        Ω = VLᵀmat(joint.qoff) * Lᵀmat(ωbar(body1, Δt)) * Lᵀmat(body1.q[No]) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end


@inline function ∂g∂posb(joint::Rotational3{T}, body1::Origin, body2::Body, No) where T
    if body2.id == joint.cid
        X = @SMatrix zeros(T, 3, 3)
        R = VLᵀmat(joint.qoff) * LVᵀmat(body2.q[No])

        return [X R]
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂velb(joint::Rotational3{T}, body1::Origin, body2::Body, Δt, No) where T
    if body2.id == joint.cid
        V = @SMatrix zeros(T, 3, 3)
        Ω = VLᵀmat(joint.qoff) * Lmat(body2.q[No]) * derivωbar(body2, Δt)

        return [V Ω]
    else
        return ∂g∂velb(joint)
    end
end
