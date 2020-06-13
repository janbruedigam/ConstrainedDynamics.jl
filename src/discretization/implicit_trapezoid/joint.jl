# Translational

@inline function g(joint::Translational, xa, qa, xb, qb)
    vertices = joint.vertices
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end
@inline function g(joint::Translational, xb, qb)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end

@inline function g(joint::Translational, statea::State, stateb::State, Δt)
    return g(joint, getx3(statea, Δt), getq3(statea, Δt), getx3(stateb, Δt), getq3(stateb, Δt))
end
@inline function g(joint::Translational, stateb::State, Δt)
    return g(joint, getx3(stateb, Δt), getq3(stateb, Δt))
end


@inline function ∂g∂posa(joint::Translational, xa, qa, xb, qb)
    point2 = xb + vrotate(joint.vertices[2], qb)

    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * (Lmat(UnitQuaternion(point2)) - Lmat(UnitQuaternion(xa))) * LVᵀmat(qa)

    return [X Q]
end
@inline function ∂g∂posb(joint::Translational, qa, qb)
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2])) * LVᵀmat(qb)

    return [X Q]
end
@inline function ∂g∂posb(joint::Translational{T}, qb) where T
    X = SMatrix{3,3,T,9}(I)
    Q = 2 * VRᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2])) * LVᵀmat(qb)

    return [X Q]
end

@inline function ∂g∂posa(joint::Translational, statea::State, stateb::State)
    return ∂g∂posa(joint, statea.xk[2], statea.qk[2], stateb.xk[2], stateb.qk[2])
end
@inline function ∂g∂posb(joint::Translational, statea::State, stateb::State)
    return ∂g∂posb(joint, statea.qk[2], stateb.qk[2])
end
@inline function ∂g∂posb(joint::Translational, stateb::State)
    return ∂g∂posb(joint, stateb.qk[2])
end


@inline function ∂g∂vela(joint::Translational, xa2, qa1, qa2, ωa, xb2, qb2, Δt)
    point2 = xb2 + vrotate(joint.vertices[2], qb2)

    V = -VLᵀmat(qa2) * RVᵀmat(qa2) * Δt
    W = 2 * VLᵀmat(qa2) * (Lmat(UnitQuaternion(point2)) - Lmat(UnitQuaternion(xa2))) * Lmat(qa1)*derivωbar(ωa, Δt)

    return [V W]
end
@inline function ∂g∂velb(joint::Translational, qa2, qb1, qb2, ωb, Δt)
    V = VLᵀmat(qa2) * RVᵀmat(qa2) * Δt
    W = 2 * VLᵀmat(qa2) * Rmat(qa2) * Rᵀmat(qb2) * Rmat(UnitQuaternion(joint.vertices[2])) * Lmat(qb1)*derivωbar(ωb, Δt)

    return [V W]
end
@inline function ∂g∂velb(joint::Translational{T}, qb1, qb2, ωb, Δt) where T
    V = SMatrix{3,3,T,9}(I) * Δt
    W = 2 * VRᵀmat(qb2) * Rmat(UnitQuaternion(joint.vertices[2])) * Lmat(qb1)*derivωbar(ωb, Δt)

    return [V W]
end

@inline function ∂g∂vela(joint::Translational, statea::State, stateb::State, Δt)
    return ∂g∂vela(joint, getx3(statea,Δt), statea.qk[2], getq3(statea,Δt), statea.ωsol[2], getx3(stateb,Δt), getq3(stateb,Δt), Δt)
end
@inline function ∂g∂velb(joint::Translational, statea::State, stateb::State, Δt)
    return ∂g∂velb(joint, getq3(statea,Δt), stateb.qk[2], getq3(stateb,Δt), stateb.ωsol[2], Δt)
end
@inline function ∂g∂velb(joint::Translational, stateb::State, Δt)
    return ∂g∂velb(joint, stateb.qk[2], getq3(stateb,Δt), stateb.ωsol[2], Δt)
end


# Rotational

@inline function g(joint::Rotational, qa, qb)
    return joint.qoff \ (qa \ qb)
end
@inline function g(joint::Rotational, qb)
    return joint.qoff \ qb
end

@inline function g(joint::Rotational, statea::State, stateb::State, Δt)
    return Vmat(g(joint, getq3(statea,Δt), getq3(stateb,Δt)))
end
@inline function g(joint::Rotational, stateb::State, Δt)
    return Vmat(g(joint, getq3(stateb,Δt)))
end


@inline function ∂g∂posa(joint::Rotational{T}, qa, qb) where T
    X = @SMatrix zeros(T, 3, 3)
    # R = VLᵀmat(joint.qoff) * Rmat(qb) * Tmat() * LVᵀmat(qa)
    R = -VLᵀmat(joint.qoff) * Rmat(qb) * RᵀVᵀmat(qa)

    return [X R]
end
@inline function ∂g∂posb(joint::Rotational{T}, qa, qb) where T
    X = @SMatrix zeros(T, 3, 3)
    R = VLᵀmat(joint.qoff) * Lᵀmat(qa) * LVᵀmat(qb)

    return [X R]
end
@inline function ∂g∂posb(joint::Rotational{T}, qb) where T
    X = @SMatrix zeros(T, 3, 3)
    R = VLᵀmat(joint.qoff) * LVᵀmat(qb)

    return [X R]
end

@inline function ∂g∂posa(joint::Rotational, statea::State, stateb::State)
    return ∂g∂posa(joint, statea.qk[2], stateb.qk[2])
end
@inline function ∂g∂posb(joint::Rotational, statea::State, stateb::State)
    return ∂g∂posb(joint, statea.qk[2], stateb.qk[2])
end
@inline function ∂g∂posb(joint::Rotational, stateb::State)
    return ∂g∂posb(joint,  stateb.qk[2])
end


@inline function ∂g∂vela(joint::Rotational{T}, qa1, ωa, qb2, Δt) where T
    V = (@SMatrix zeros(T, 3, 3)) * Δt
    W = VLᵀmat(joint.qoff) * Rmat(qb2) * Tmat(T) * Lmat(qa1)*derivωbar(ωa, Δt)

    return [V W]
end
@inline function ∂g∂velb(joint::Rotational{T}, qa2, qb1, ωb, Δt) where T
    V = (@SMatrix zeros(T, 3, 3)) * Δt
    W = VLᵀmat(joint.qoff) * Lᵀmat(qa2) * Lmat(qb1)*derivωbar(ωb, Δt)

    return [V W]
end
@inline function ∂g∂velb(joint::Rotational{T}, qb1, ωb, Δt) where T
    V = (@SMatrix zeros(T, 3, 3)) * Δt
    W = VLᵀmat(joint.qoff) * Lmat(qb1)*derivωbar(ωb, Δt)

    return [V W]
end

@inline function ∂g∂vela(joint::Rotational, statea::State, stateb::State, Δt)
    return ∂g∂vela(joint, statea.qk[2], statea.ωsol[2], getq3(stateb,Δt), Δt)
end
@inline function ∂g∂velb(joint::Rotational, statea::State, stateb::State, Δt)
    return ∂g∂velb(joint, getq3(statea,Δt), stateb.qk[2], stateb.ωsol[2], Δt)
end
@inline function ∂g∂velb(joint::Rotational, stateb::State, Δt)
    return ∂g∂velb(joint, stateb.qk[2], stateb.ωsol[2], Δt)
end


# Forcing

@inline function setForce!(joint::Translational, statea::State, stateb::State, F)
    vertices = joint.vertices
    qa = statea.qk[2]
    qb = stateb.qk[2]

    Fa = vrotate(-F, qa)
    Fb = -Fa

    τa = vrotate(torqueFromForce(Fa, vrotate(vertices[1], qa)),inv(qa)) # in local coordinates
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    statea.Fk[2] += Fa
    statea.τk[2] += τa
    stateb.Fk[2] += Fb
    stateb.τk[2] += τb
    return
end
@inline function setForce!(joint::Translational, stateb::State, F)
    vertices = joint.vertices
    qb = stateb.qk[2]

    Fb = F
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    stateb.Fk[2] += Fb
    stateb.τk[2] += τb
    return
end

@inline function setForce!(joint::Rotational, statea::State, stateb::State, τ)
    qa = statea.qk[2]
    qb = stateb.qk[2]    

    τa = vrotate(-τ, qa*joint.qoff) # in world coordinates
    τb = -τa # in world coordinates

    τa = vrotate(τa,inv(qa)) # in local coordinates
    τb = vrotate(τb,inv(qb)) # in local coordinates

    statea.τk[2] += τa
    stateb.τk[2] += τb
    return
end
@inline function setForce!(joint::Rotational, stateb::State, τ)
    qb = stateb.qk[2]

    τa = vrotate(-τ, joint.qoff) # in world coordinates
    τb = -τa # in world coordinates

    τb = vrotate(τb,inv(qb)) # in local coordinates

    stateb.τk[2] += τb
    return
end


@inline function ∂Fτ∂ua(joint::Translational, statea::State)
    vertices = joint.vertices
    qa = statea.qk[2]

    BFa = -VLmat(qa) * RᵀVᵀmat(qa)
    Bτa = VLᵀmat(qa) * RVᵀmat(qa) * skew(BFa*vertices[1])

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Translational, statea::State, stateb::State)
    vertices = joint.vertices
    qa = statea.qk[2]
    qb = stateb.qk[2]

    BFb = VLmat(qa) * RᵀVᵀmat(qa)
    Bτb = VLᵀmat(qb) * RVᵀmat(qb) * skew(BFb*vertices[2])

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Translational{T}, stateb::State) where T
    vertices = joint.vertices
    qb = stateb.qk[2]

    BFb = SMatrix{3,3,T,9}(I)
    Bτb = VLᵀmat(qb) * RVᵀmat(qb) * skew(vertices[2])

    return [BFb; Bτb]
end

@inline function ∂Fτ∂ua(joint::Rotational{T}, statea::State) where T
    qoff = joint.qoff

    BFa = (@SMatrix zeros(T, 3, 3))
    Bτa = -VLmat(qoff) * RᵀVᵀmat(qoff)

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, statea::State, stateb::State) where T
    qa = statea.qk[2]
    qb = stateb.qk[2]
    qoff = joint.qoff
    qaiqaqoff = qb\qa*qoff

    BFb = (@SMatrix zeros(T, 3, 3))
    Bτb = VLmat(qaiqaqoff) * RᵀVᵀmat(qaiqaqoff)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, stateb::State) where T
    qb = stateb.qk[2]
    qoff = joint.qoff
    qaiqoff = qb\qoff

    BFb = (@SMatrix zeros(T, 3, 3))
    Bτb = VLmat(qaiqoff) * RᵀVᵀmat(qaiqoff)

    return [BFb; Bτb]
end