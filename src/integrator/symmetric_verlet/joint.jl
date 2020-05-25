# Translational

@inline function g(joint::Translational, xa, qa, xb, qb)
    vertices = joint.vertices
    vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end

@inline function g(joint::Translational, xb, qb)
    vertices = joint.vertices
    xb + vrotate(vertices[2], qb) - vertices[1]
end


@inline function ∂g∂posa(joint::Translational, xa, qa, xb, qb)
    point2 = xb + vrotate(joint.vertices[2], qb)

    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(xa))) * LVᵀmat(qa)

    return [X Q]
end

@inline function ∂g∂posb(joint::Translational, xa, qa, xb, qb)
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(qb)

    return [X Q]
end

@inline function ∂g∂posb(joint::Translational{T}, xb, qb) where T
    X = SMatrix{3,3,T,9}(I)
    Q = 2 * VRᵀmat(qb) * Rmat(Quaternion(joint.vertices[2])) * LVᵀmat(qb)

    return [X Q]
end


@inline function ∂g∂vela(joint::Translational, xa1, xa2, qa1, qa2, va2, ωa2, xb2, qb2, Δt)
    point2 = xb2 + vrotate(joint.vertices[2], qb2)

    V = -VLᵀmat(qa2) * RVᵀmat(qa2) * 1/2 * Δt
    W = 2 * VLᵀmat(qa2) * (Lmat(Quaternion(point2)) - Lmat(Quaternion(xa2))) * 1/2 * Lmat(qa1)*derivωbar(ωa2, Δt)

    return [V W]
end

@inline function ∂g∂velb(joint::Translational, xa2, qa2, xb1, xb2, qb1, qb2, vb2, ωb2, Δt)
    V = VLᵀmat(qa2) * RVᵀmat(qa2) * 1/2 * Δt
    W = 2 * VLᵀmat(qa2) * Rmat(qa2) * Rᵀmat(qb2) * Rmat(Quaternion(joint.vertices[2])) * 1/2 * Lmat(qb1)*derivωbar(ωb2, Δt)

    return [V W]
end


@inline function ∂g∂velb(joint::Translational{T}, xb1, xb2, qb1, qb2, vb2, ωb2, Δt) where T
    V = SMatrix{3,3,T,9}(I) * 1/2 * Δt
    W = 2 * VRᵀmat(qb2) * Rmat(Quaternion(joint.vertices[2])) * 1/2 * Lmat(qb1)*derivωbar(ωb2, Δt)

    return [V W]
end



# Rotational

@inline function g(joint::Rotational, xa, qa, xb, qb)
    VLᵀmat(joint.qoff) * Lᵀmat(qa) * qb
end

@inline function g(joint::Rotational, xb, qb)
    VLᵀmat(joint.qoff) * qb
end


@inline function ∂g∂posa(joint::Rotational{T}, xa, qa, xb, qb) where T
    X = @SMatrix zeros(T, 3, 3)
    # R = VLᵀmat(joint.qoff) * Rmat(qb) * Tmat() * LVᵀmat(qa)
    R = -VLᵀmat(joint.qoff) * Rmat(qb) * RᵀVᵀmat(qa)

    return [X R]
end

@inline function ∂g∂posb(joint::Rotational{T}, xa, qa, xb, qb) where T
    X = @SMatrix zeros(T, 3, 3)
    R = VLᵀmat(joint.qoff) * Lᵀmat(qa) * LVᵀmat(qb)

    return [X R]
end

@inline function ∂g∂posb(joint::Rotational{T}, xb, qb) where T
    X = @SMatrix zeros(T, 3, 3)
    R = VLᵀmat(joint.qoff) * LVᵀmat(qb)

    return [X R]
end


@inline function ∂g∂vela(joint::Rotational{T}, xa1, xa2, qa1, qa2, va2, ωa2, xb2, qb2, Δt) where T
    V = (@SMatrix zeros(T, 3, 3)) * 1/2 * Δt
    W = VLᵀmat(joint.qoff) * Rmat(qb2) * Tmat(T) * 1/2 * Lmat(qa1)*derivωbar(ωa2, Δt)

    return [V W]
end

@inline function ∂g∂velb(joint::Rotational{T}, xa2, qa2, xb1, xb2, qb1, qb2, vb2, ωb2, Δt) where T
    V = (@SMatrix zeros(T, 3, 3)) * 1/2 * Δt
    W = VLᵀmat(joint.qoff) * Lᵀmat(qa2) * 1/2 * Lmat(qb1)*derivωbar(ωb2, Δt)

    return [V W]
end


@inline function ∂g∂velb(joint::Rotational{T}, xb1, xb2, qb1, qb2, vb2, ωb2, Δt) where T
    V = (@SMatrix zeros(T, 3, 3)) * 1/2 * Δt
    W = VLᵀmat(joint.qoff) * 1/2 * Lmat(qb1)*derivωbar(ωb2, Δt)

    return [V W]
end