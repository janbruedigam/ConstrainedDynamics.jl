# Forcing

@inline function setForce!(joint::Translational, statea::State, stateb::State, F)
    vertices = joint.vertices
    qa = statea.qk[1]
    qb = stateb.qk[1]

    Fa = vrotate(-F, qa)
    Fb = -Fa

    τa = vrotate(torqueFromForce(Fa, vrotate(vertices[1], qa)),inv(qa)) # in local coordinates
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    statea.Fk[1] += Fa
    statea.τk[1] += τa
    stateb.Fk[1] += Fb
    stateb.τk[1] += τb
    return
end
@inline function setForce!(joint::Translational, stateb::State, F)
    vertices = joint.vertices
    qb = stateb.qk[1]

    Fb = F
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    stateb.Fk[1] += Fb
    stateb.τk[1] += τb
    return
end

@inline function setForce!(joint::Rotational, statea::State, stateb::State, τ)
    qa = statea.qk[1]
    qb = stateb.qk[1]    

    τa = vrotate(-τ, qa*joint.qoff) # in world coordinates
    τb = -τa # in world coordinates

    τa = vrotate(τa,inv(qa)) # in local coordinates
    τb = vrotate(τb,inv(qb)) # in local coordinates

    statea.τk[1] += τa
    stateb.τk[1] += τb
    return
end
@inline function setForce!(joint::Rotational, stateb::State, τ)
    qb = stateb.qk[1]

    τa = vrotate(-τ, joint.qoff) # in world coordinates
    τb = -τa # in world coordinates

    τb = vrotate(τb,inv(qb)) # in local coordinates

    stateb.τk[1] += τb
    return
end


# Control derivatives

@inline function ∂Fτ∂ua(joint::Translational, statea::State)
    vertices = joint.vertices
    qa = statea.qk[1]

    BFa = -VLmat(qa) * RᵀVᵀmat(qa)
    Bτa = VLᵀmat(qa) * RVᵀmat(qa) * skew(BFa*vertices[1])

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Translational, statea::State, stateb::State)
    vertices = joint.vertices
    qa = statea.qk[1]
    qb = stateb.qk[1]

    BFb = VLmat(qa) * RᵀVᵀmat(qa)
    Bτb = VLᵀmat(qb) * RVᵀmat(qb) * skew(BFb*vertices[2])

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Translational{T}, stateb::State) where T
    vertices = joint.vertices
    qb = stateb.qk[1]

    BFb = I
    Bτb = VLᵀmat(qb) * RVᵀmat(qb) * skew(vertices[2])

    return [BFb; Bτb]
end

@inline function ∂Fτ∂ua(joint::Rotational{T}, statea::State) where T
    qoff = joint.qoff

    BFa = (szeros(T, 3, 3))
    Bτa = -VLmat(qoff) * RᵀVᵀmat(qoff)

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, statea::State, stateb::State) where T
    qa = statea.qk[1]
    qb = stateb.qk[1]
    qoff = joint.qoff
    qaiqaqoff = qb\qa*qoff

    BFb = (szeros(T, 3, 3))
    Bτb = VLmat(qaiqaqoff) * RᵀVᵀmat(qaiqaqoff)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, stateb::State) where T
    qb = stateb.qk[1]
    qoff = joint.qoff
    qaiqoff = qb\qoff

    BFb = (szeros(T, 3, 3))
    Bτb = VLmat(qaiqoff) * RᵀVᵀmat(qaiqoff)

    return [BFb; Bτb]
end