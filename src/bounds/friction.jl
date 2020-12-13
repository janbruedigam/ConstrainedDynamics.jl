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


@inline function frictionmap(friction::Friction, x::AbstractVector, q::UnitQuaternion)
    return [I;VLᵀmat(q)*RVᵀmat(q)*skew(friction.p)] * friction.D'
end
@inline function invfrictionmap(friction::Friction, x::AbstractVector, q::UnitQuaternion)
    return friction.D * [I [friction.p';friction.p';friction.p']*VRᵀmat(q)*LVᵀmat(q)]
end

@inline additionalforce(friction::Friction, x::AbstractVector, q::UnitQuaternion) = frictionmap(friction, x, q)*friction.b

@inline function calcFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    calcFrictionForce!(mechanism, friction, body, ineqc.γsol[2][i])
    return
end

##Friction
@inline function calcFrictionForce!(mechanism::Mechanism{T}, friction::Contact, body::Body, γ) where T
    cf = friction.cf
    D = friction.D
    state = body.state

    d = state.d
    v = state.vsol[2]
    ω = state.ωsol[2]
    state.vsol[2] = szeros(T, 3)
    state.ωsol[2] = szeros(T, 3)
    dyn = dynamics(mechanism, body)
    state.vsol[2] = v
    state.ωsol[2] = ω
    state.d = d

    b0 = invfrictionmap(friction, posargsk(body.state)...)*dyn + friction.b # remove old friction force

    if norm(b0) > cf*γ
        friction.b = b0/norm(b0)*cf*γ
    else
        friction.b = b0
    end    
    return
end