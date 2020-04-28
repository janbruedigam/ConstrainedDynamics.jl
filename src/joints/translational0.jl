@inline function getVelocityDelta(joint::Translational0, body1::AbstractBody, body2::Body{T}, v::SVector{3,T}) where T
    Δv = v # in body1 frame
    return Δv
end

@inline function getPositionDelta(joint::Translational0, body1::AbstractBody, body2::Body{T}, x::SVector{3,T}) where T
    Δx = x # in body1 frame
    return Δx
end

@inline function setForce!(joint::Translational0, body1::Body, body2::Body{T}, F::SVector{3,T}, No) where T
    F1 = vrotate(-F, body1.q[No])
    F2 = -F1

    body1.F[No] = F1
    body2.F[No] = F2
    return
end

@inline function setForce!(joint::Translational0, body1::Origin, body2::Body{T}, F::SVector{3,T}, No) where T
    body2.F[No] = F
    return
end


@inline function minimalCoordinates(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No)
    body2.x[No]
end


@inline g(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = g(joint)


@inline ∂g∂posa(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Translational0, body1::AbstractBody, body2::AbstractBody, No) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Translational0, body1::AbstractBody, body2::AbstractBody, Δt, No) = ∂g∂velb(joint)
