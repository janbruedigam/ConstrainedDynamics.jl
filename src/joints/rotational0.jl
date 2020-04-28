@inline function getVelocityDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, ω::SVector{3,T}) where T
    Δω = ω # in body1 frame
    return Δω
end

@inline function getPositionDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, θ::Quaternion{T}) where T
    Δq = θ # in body1 frame
    return Δq
end

@inline function setForce!(joint::Rotational0, body1::Body, body2::Body{T}, τ::SVector{3,T}, No) where T
    τ1 = vrotate(-τ, body1.q[No] * joint.qoff)
    τ2 = -τ1

    body1.τ[No] = τ1
    body2.τ[No] = τ2
    return
end

@inline function setForce!(joint::Rotational0, body1::Origin, body2::Body{T}, τ::SVector{3,T}, No) where T
    body2.τ[No] = vrotate(τ, joint.qoff)
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
