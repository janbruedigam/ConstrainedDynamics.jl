abstract type Joint{T,N} <: AbstractJoint{T,N} end

Joint0 = Joint{T,0} where T
Joint1 = Joint{T,1} where T
Joint2 = Joint{T,2} where T
Joint3 = Joint{T,3} where T

Base.show(io::IO, joint::Joint) = summary(io, joint)

getT(joint::AbstractJoint{T}) where T = T
getN(joint::AbstractJoint{T,N}) where {T,N} = N

@inline constraintmat(::Joint0{T}) where T = szeros(T,0,3)
@inline nullspacemat(::Joint0{T}) where T = SMatrix{3,3,T,9}(I)
@inline constraintmat(joint::Joint1) = joint.V3
@inline nullspacemat(joint::Joint1) = joint.V12
@inline constraintmat(joint::Joint2) = joint.V12
@inline nullspacemat(joint::Joint2) = joint.V3
@inline constraintmat(::Joint3{T}) where T = SMatrix{3,3,T,9}(I)
@inline nullspacemat(::Joint3{T}) where T = szeros(T,0,3)

@inline setForce!(joint::Joint) = return

@inline minimalCoordinates(joint::Joint{T,N}) where {T,N} = szeros(T, 3 - N)
@inline g(joint::Joint{T,N}) where {T,N} = szeros(T, N)

@inline ∂g∂ʳposa(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ʳposb(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ʳvela(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 6)
@inline ∂g∂ʳvelb(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 6)

@inline ∂g∂posac(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 7)
@inline ∂g∂posbc(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 7)

@inline ∂Fτ∂ub(joint::Joint{T,N}) where {T,N} = szeros(T, 6, 3 - N)

@inline function setForce!(joint::Joint, body1::Body, body2::Body, Fτ::SVector)
    setForce!(joint, body1.state, body2.state, zerodimstaticadjoint(nullspacemat(joint)) * Fτ)
    return
end
@inline function setForce!(joint::Joint, body1::Origin, body2::Body, Fτ::SVector)
    setForce!(joint, body2.state, zerodimstaticadjoint(nullspacemat(joint)) * Fτ)
    return
end

@inline g(joint::AbstractJoint, body1::Body, body2::Body, Δt) = constraintmat(joint) * g(joint, body1.state, body2.state, Δt)
@inline g(joint::AbstractJoint, body1::Origin, body2::Body, Δt) = constraintmat(joint) * g(joint, body2.state, Δt)
@inline g(joint::AbstractJoint, body1::Body, body2::Body) = constraintmat(joint) * g(joint, body1.state, body2.state)
@inline g(joint::AbstractJoint, body1::Origin, body2::Body) = constraintmat(joint) * g(joint, body2.state)

@inline function ∂g∂ʳposa(joint::AbstractJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂ʳposa(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposa(joint)
    end
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂ʳposb(joint, body1.state, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂ʳposb(joint, body2.state)
    else
        return ∂g∂ʳposb(joint)
    end
end

@inline function ∂g∂posac(joint::AbstractJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂posac(joint, body1.state, body2.state)
    else
        return ∂g∂posac(joint)
    end
end
@inline function ∂g∂posbc(joint::AbstractJoint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂posbc(joint, body1.state, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end
@inline function ∂g∂posbc(joint::AbstractJoint, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂posbc(joint, body2.state)
    else
        return ∂g∂posbc(joint)
    end
end

@inline function ∂g∂ʳvela(joint::AbstractJoint, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvela(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, body1::Body, body2::Body, Δt)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, body1::Origin, body2::Body, Δt)
    if body2.id == joint.childid
        return constraintmat(joint) * ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return ∂g∂ʳvelb(joint)
    end
end

@inline function ∂Fτ∂ua(joint::Joint, body1::Body)
    return ∂Fτ∂ua(joint, body1.state) * zerodimstaticadjoint(nullspacemat(joint))
end
@inline function ∂Fτ∂ub(joint::Joint, body1::Body, body2::Body)
    if body2.id == joint.childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Joint, body1::Origin, body2::Body)
    if body2.id == joint.childid
        return return ∂Fτ∂ub(joint, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end
end


# Derivatives accounting for quaternion specialness

@inline function ∂g∂ʳposa(joint::AbstractJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posa(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qa)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xa, qa, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, xb::AbstractVector, qb::UnitQuaternion)
    X, Q = ∂g∂posb(joint, xb, qb)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end

@inline function ∂g∂ʳvela(joint::AbstractJoint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, Δt
    )

    X, Q = ∂g∂posa(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1a) * derivωbar(ω1a, Δt)

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
    x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2a, q2a, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt)

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, Δt
    )

    X, Q = ∂g∂posb(joint, x2b, q2b)
    V = X * Δt
    Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt)

    return [V Ω]
end


# Calls for dynamics
g(joint::Joint, statea::State, stateb::State, Δt) = g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::Joint, stateb::State, Δt) = g(joint, posargsnext(stateb, Δt)...)
g(joint::Joint, statea::State, stateb::State) = g(joint, posargsc(statea)..., posargsc(stateb)...)
g(joint::Joint, stateb::State) = g(joint, posargsc(stateb)...)

∂g∂ʳposa(joint::AbstractJoint, statea::State, stateb::State) = ∂g∂ʳposa(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::AbstractJoint, statea::State, stateb::State) = ∂g∂ʳposb(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::AbstractJoint, stateb::State) = ∂g∂ʳposb(joint, posargsk(stateb)...)

∂g∂ʳvela(joint::AbstractJoint, statea::State, stateb::State, Δt) = ∂g∂ʳvela(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(statea)..., Δt)
∂g∂ʳvelb(joint::AbstractJoint, statea::State, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)
∂g∂ʳvelb(joint::AbstractJoint, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)

# Calls for constraints solver
∂g∂posac(joint::AbstractJoint, statea::State, stateb::State) = hcat(∂g∂posa(joint, posargsc(statea)..., posargsc(stateb)...)...)
∂g∂posbc(joint::AbstractJoint, statea::State, stateb::State) = hcat(∂g∂posb(joint, posargsc(statea)..., posargsc(stateb)...)...)
∂g∂posbc(joint::AbstractJoint, stateb::State) = hcat(∂g∂posb(joint, posargsc(stateb)...)...)
