mutable struct Friction{T} <: Contact{T}
    ainv3::Adjoint{T,SVector{3,T}}
    p::SVector{3,T}
    D::SMatrix{2,3,T,6}
    cf::T
    offset::SVector{3,T}
    b::SVector{2,T}
    

    function Friction(body::Body{T}, normal::AbstractVector, cf::Real; p::AbstractVector = zeros(3), offset::AbstractVector = zeros(3)) where T
        @assert cf>0

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        D = [V1 V2]'

        new{T}(ainv3, p, D, cf, offset, zeros(2)), body.id
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Friction{T}) where {T}
    summary(io, constraint)
    println(io,"")
    println(io, " ainv3:  "*string(constraint.ainv3))
    println(io, " D:      "*string(constraint.D))
    println(io, " cf:     "*string(constraint.cf))
    println(io, " offset: "*string(constraint.offset))
end

@inline function calcb0(friction::Friction, m, Δt, x::AbstractVector, q::UnitQuaternion, vc, ωc, dyn)
    p = friction.p
    np2 = norm(p)^2
    vp = vc + vrotate(cross(ωc,p),q)
    if np2 > 0
        return friction.D * (I*m*(-vp)/Δt + dyn[SA[1;2;3]] + vrotate(cross(dyn[SA[4;5;6]]/2,p),q)/np2)
    else
        return friction.D * (I*m*(-vp)/Δt + dyn[SA[1;2;3]])
    end
end

@inline function frictionmap(friction::Friction, x::AbstractVector, q::UnitQuaternion)
    Fp = friction.D' * friction.b
    Fc = Fp
    τc = torqueFromForce(vrotate(Fp,inv(q)), friction.p) # in local coordinates
    return [Fc;2*τc]
end

@inline additionalforce(friction::Friction, x::AbstractVector, q::UnitQuaternion) = frictionmap(friction, x, q)

@inline function calcFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    calcFrictionForce!(mechanism, friction, body, ineqc.γsol[2][i])
    return
end

##Friction
# TODO this might be wrong, will be done differently in future anyways
@inline function calcFrictionForce!(mechanism::Mechanism{T}, friction::Contact, body::Body, γ) where T
    cf = friction.cf
    state = body.state

    d = state.d
    v = state.vsol[2]
    ω = state.ωsol[2]
    state.vsol[2] = state.vc
    state.ωsol[2] = state.ωc
    dyn = dynamics(mechanism, body) + frictionmap(friction, posargsk(body.state)...)
    state.vsol[2] = v
    state.ωsol[2] = ω
    state.d = d

    b0 = calcb0(friction, body.m, mechanism.Δt, posargsk(body.state)..., state.vc, state.ωc, dyn)
    
    if norm(b0) > cf*γ
        friction.b = b0/norm(b0)*cf*γ
    else
        friction.b = b0
    end    
    return
end