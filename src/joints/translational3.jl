@inline function getPositionDelta(joint::Translational3, body1::AbstractBody, body2::Body, x::SVector{0,T}) where T
    Δx = @SVector zeros(T,3)
    return Δx
end
@inline function getVelocityDelta(joint::Translational3, body1::AbstractBody, body2::Body, v::SVector{0,T}) where T
    Δv = @SVector zeros(T,3)
    return Δv
end

@inline function ∂Fτ∂ua(joint::Translational3{T}, body1::Body) where T
    return @SMatrix zeros(T, 6, 0)
end
@inline function ∂Fτ∂ub(joint::Translational3{T}, body1::AbstractBody, body2::Body) where T
    return @SMatrix zeros(T, 6, 0)
end

@inline minimalCoordinates(joint::Translational3{T}, body1::AbstractBody, body2::Body) where T = SVector{0,T}()

@inline g(joint::Translational3, body1::Body, body2::Body, Δt) = g(joint, body1.state, body2.state, Δt)
@inline g(joint::Translational3, body1::Origin, body2::Body, Δt) = g(joint, body2.state, Δt)

@inline function ∂g∂posa(joint::Translational3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posa(joint, body1.state, body2.state)
    else
        return ∂g∂posa(joint)
    end
end
@inline function ∂g∂posb(joint::Translational3, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posb(joint, body1.state, body2.state)
    else
        return ∂g∂posb(joint)
    end
end
@inline function ∂g∂posb(joint::Translational3, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return ∂g∂posb(joint, body2.state)
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂vela(joint::Translational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂vela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂vela(joint)
    end
end
@inline function ∂g∂velb(joint::Translational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂velb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end
@inline function ∂g∂velb(joint::Translational3, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return ∂g∂velb(joint, body2.state, Δt)
    else
        return ∂g∂velb(joint)
    end
end