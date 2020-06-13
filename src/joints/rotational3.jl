@inline function getPositionDelta(joint::Rotational3, body1::AbstractBody, body2::Body, θ::SVector{0,T}) where T
    Δq = joint.qoff
    return Δq
end
@inline function getVelocityDelta(joint::Rotational3, body1::AbstractBody, body2::Body, ω::SVector{0,T}) where T
    Δω = @SVector zeros(T,3)
    return Δω
end

@inline function ∂Fτ∂ua(joint::Rotational3{T}, body1::Body) where T
    return @SMatrix zeros(T, 6, 0)
end
@inline function ∂Fτ∂ub(joint::Rotational3{T}, body1::AbstractBody, body2::Body) where T
    return @SMatrix zeros(T, 6, 0)
end

@inline minimalCoordinates(joint::Rotational3{T}, body1::AbstractBody, body2::Body) where T = SVector{0,T}()

@inline g(joint::Rotational3, body1::Body, body2::Body, Δt) = g(joint, body1.state, body2.state, Δt)
@inline g(joint::Rotational3, body1::Origin, body2::Body, Δt) = g(joint, body2.state, Δt)

@inline function ∂g∂posa(joint::Rotational3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posa(joint, body1.state, body2.state)
    else
        return ∂g∂posa(joint)
    end
end
@inline function ∂g∂posb(joint::Rotational3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posb(joint, body1.state, body2.state)
    else
        return ∂g∂posb(joint)
    end
end
@inline function ∂g∂posb(joint::Rotational3, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posb(joint, body2.state)
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Rotational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂vela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂vela(joint)
    end
end
@inline function ∂g∂velb(joint::Rotational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂velb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
@inline function ∂g∂velb(joint::Rotational3, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂velb(joint, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end