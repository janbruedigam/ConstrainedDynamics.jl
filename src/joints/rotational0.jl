@inline function getVelocityDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, ω::SVector{3,T}) where T
    Δω = ω # in body1 frame
    return Δω
end

@inline function getPositionDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, θ::Quaternion{T}) where T
    Δq = θ # in body1 frame
    return Δq
end

@inline function setForce!(joint::Rotational0, body1::Body, body2::Body{T}, τ::SVector{3,T}, No) where T
    clearForce!(joint, body1, body2, No)

    τ1 = vrotate(-τ, body1.q[No])
    τ2 = -τ1
    F1 =  @SVector zeros(T,3)
    F2 =  @SVector zeros(T,3)

    updateForce!(joint, body1, body2, F1, τ1, F2, τ2, No)
    return
end

@inline function setForce!(joint::Rotational0, body1::Origin, body2::Body{T}, τ::SVector{3,T}, No) where T
    clearForce!(joint, body2, No)

    τ2 = τ
    F2 =  @SVector zeros(T,3)

    updateForce!(joint, body2, F2, τ2, No)
    return
end


@inline function minimalCoordinates(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No)
    body2.q[No]
end


@inline g(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = g(joint)


@inline ∂g∂posa(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂velb(joint)
