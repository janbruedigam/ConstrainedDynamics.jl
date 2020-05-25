@inline function getPositionDelta(joint::Translational3, body1::AbstractBody, body2::Body{T}, x::SVector{0,T}) where T
    Δx = @SVector zeros(T,3)
    return Δx
end

@inline function getVelocityDelta(joint::Translational3, body1::AbstractBody, body2::Body{T}, v::SVector{0,T}) where T
    Δv = @SVector zeros(T,3)
    return Δv
end

@inline function minimalCoordinates(joint::Translational3, body1::AbstractBody{T}, body2::Body, No) where T
    SVector{0,T}()
end


@inline function g(joint::Translational3, body1::Body, body2::Body, Δt)
    g(joint, getx2(body1, Δt), getq2(body1, Δt), getx2(body2, Δt), getq2(body2, Δt))
end

@inline function g(joint::Translational3, body1::Origin, body2::Body, Δt)
    g(joint, getx2(body2, Δt), getq2(body2, Δt))
end


@inline function ∂g∂posa(joint::Translational3, body1::Body, body2::Body)
    if body2.id == joint.cid
        return ∂g∂posa(joint, getxd2(body1), getqd2(body1), getxd2(body2), getqd2(body2))
    else
        return ∂g∂posa(joint)
    end
end

@inline function ∂g∂posb(joint::Translational3, body1::Body, body2::Body)
    if body2.id == joint.cid
        return ∂g∂posb(joint, getxd2(body1), getqd2(body1), getxd2(body2), getqd2(body2))
    else
        return ∂g∂posb(joint)
    end
end

@inline function ∂g∂posb(joint::Translational3, body1::Origin, body2::Body)
    if body2.id == joint.cid
        return ∂g∂posb(joint, getxd2(body2), getqd2(body2))
    else
        return ∂g∂posb(joint)
    end
end


@inline function ∂g∂vela(joint::Translational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return ∂g∂vela(joint, getxd2(body1), getx2(body1, Δt), getqd2(body1), getq2(body1, Δt), getv2(body1), getω2(body1), getx2(body2, Δt), getq2(body2, Δt), Δt)
    else
        return ∂g∂vela(joint)
    end
end

@inline function ∂g∂velb(joint::Translational3, body1::Body, body2::Body, Δt)
    if body2.id == joint.cid
        return ∂g∂velb(joint, getx2(body1, Δt), getq2(body1, Δt), getxd2(body2), getx2(body2, Δt), getqd2(body2), getq2(body2, Δt), getv2(body2), getω2(body2), Δt)
    else
        return ∂g∂velb(joint)
    end
end

@inline function ∂g∂velb(joint::Translational3, body1::Origin, body2::Body, Δt)
    if body2.id == joint.cid
        return ∂g∂velb(joint, getxd2(body2), getx2(body2, Δt), getqd2(body2), getq2(body2, Δt), getv2(body2), getω2(body2), Δt)
    else
        return ∂g∂velb(joint)
    end
end