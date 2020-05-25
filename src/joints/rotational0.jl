@inline function getPositionDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, θ::SVector{3,T}) where T
    #TODO define this function (choose 3d parameters)
    @error("not defined for rot0")
end

@inline function getVelocityDelta(joint::Rotational0, body1::AbstractBody, body2::Body{T}, ω::SVector{3,T}) where T
    Δω = vrotate(ω, inv(body2.state.qd[2])*body1.state.qd[2]*joint.qoff) # in body2 frame
    return Δω
end

@inline function setForce!(joint::Rotational0, body1::Body, body2::Body{T}, τ::SVector{3,T}, No) where T
    clearForce!(joint, body1, body2, No)

    q1 = body1.state.qd[No]
    q2 = body2.state.qd[No]    

    τ1 = vrotate(-τ, q1*joint.qoff) # in world coordinates
    τ2 = -τ1 # in world coordinates

    τ1 = vrotate(τ1,inv(q1)) # in local coordinates
    τ2 = vrotate(τ2,inv(q2)) # in local coordinates
    
    F1 =  @SVector zeros(T,3)
    F2 =  @SVector zeros(T,3)

    updateForce!(joint, body1, body2, F1, τ1, F2, τ2, No)
    return
end

@inline function setForce!(joint::Rotational0, body1::Origin, body2::Body{T}, τ::SVector{3,T}, No) where T
    clearForce!(joint, body2, No)

    q2 = body2.state.qd[No]

    τ1 = vrotate(-τ, joint.qoff) # in world coordinates
    τ2 = -τ1 # in world coordinates

    τ2 = vrotate(τ2,inv(q2)) # in local coordinates

    F2 =  @SVector zeros(T,3)

    updateForce!(joint, body2, F2, τ2, No)
    return
end


@inline function minimalCoordinates(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, No)
    body2.state.qd[No]
end


@inline g(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = g(joint)


@inline ∂g∂posa(joint::Rotational0, body1::AbstractBody, body2::AbstractBody) = ∂g∂posa(joint)
@inline ∂g∂posb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody) = ∂g∂posb(joint)
@inline ∂g∂vela(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂vela(joint)
@inline ∂g∂velb(joint::Rotational0, body1::AbstractBody, body2::AbstractBody, Δt) = ∂g∂velb(joint)
