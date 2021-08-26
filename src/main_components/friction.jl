mutable struct Friction{T,N} <: Component{T}
    id::Int64
    name::String

    Bx::SMatrix{4,3,T,12}
    p::SVector{3,T}
    cf::T

    parentid::Union{Int64,Nothing}
    childids::SVector{2,Int64}

    βsol::Vector{SVector{N,T}}

    function Friction(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3), name::String="") where T
        
        Bx = SA{T}[
            1 0 0
            -1 0 0
            0 1 0
            0 -1 0
        ]

        frictionid = getGlobalID()
        impact = InequalityConstraint(Impact(body, normal; p = p, offset = offset))
        frictionbound = InequalityConstraint(FrictionBound(frictionid, impact.id))
        betabound = InequalityConstraint(BetaBound(frictionid))

        new{T,4}(frictionid, name, Bx, p, cf, body.id, [frictionbound.id; betabound.id]), [impact; frictionbound; betabound]
    end
end


@inline ∂g∂ʳvel(mechanism, fric::Friction{T,N}) where {T,N} = szeros(T,N,N)

# @inline Bq(Bxmat, p, q) = Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p)

# ### Constraints and derivatives
# ## Position level constraint wrappers
# g(friction::Friction, body::Body, ineqc1, ineqc2) = g(friction, body.state, ineqc1, ineqc2)

# ## Discrete-time position wrappers (for dynamics)
# @inline g(friction::Friction, state::State, ineqc1, ineqc2) = g(friction, fullargssol(state)..., ineqc1.γsol[2], ineqc2.γsol[2])

# ## Position level constraints (for dynamics)
# @inline function g(friction::Friction, x::AbstractVector, v::AbstractVector, q::UnitQuaternion, ω::AbstractVector, ψ, η)
#     Bxmat = friction.Bx
#     Bqmat = Bq(Bxmat, friction.p, q)
#     Bxmat*v + Bqmat*ω - η .+ ψ
# end

# ## Discrete-time position derivatives (for dynamics)
# # Wrappers 1
# @inline function ∂g∂ʳpos(friction::Friction, body::Body)
#     return ∂g∂ʳpos(friction, body.state)
# end
# @inline function ∂g∂ʳpos(friction::Friction, ineqc::InequalityConstraint)
#     return ∂g∂ʳpos(friction, ineqc.constraints[1])
# end

# # Wrappers 2
# ∂g∂ʳpos(friction::Friction, state::State) = ∂g∂ʳpos(friction, posargsk(state)...)
# ∂g∂ʳpos(friction::Friction{T,N}, constraint::FrictionBound) where {T,N} = [szeros(T,N) sones(T,N)]
# ∂g∂ʳpos(friction::Friction{T,N}, constraint::BetaBound) where {T,N} = [szeros(T,N,N) -I]

# # Derivatives accounting for quaternion specialness
# @inline function ∂g∂ʳpos(friction::Friction, x::AbstractVector, q::UnitQuaternion)
#     Bxmat = friction.Bx
#     Bqmat = Bq(Bxmat, friction.p, q)

#     return [Bxmat Bqmat]
# end






# function g(mechanism, eqc::EqualityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{<:Friction}}
#     return g(eqc.constraints[1], getbody(mechanism, eqc.parentid), getineqconstraint(mechanism, eqc.childids[1]), getineqconstraint(mechanism, eqc.childids[2]))
# end
# @inline function ∂g∂ʳvel(mechanism, eqc::EqualityConstraint{T,N,Nc,Cs}, id::Integer) where {T,N,Nc,Cs<:Tuple{<:Friction}}
#     return ∂g∂ʳpos(mechanism, eqc, id)
# end
# function ∂g∂ʳposa(mechanism, eqc::EqualityConstraint{T,N,Nc,Cs}, id::Integer) where {T,N,Nc,Cs<:Tuple{<:Friction}}
#     return ∂g∂ʳpos(eqc.constraints[1], getbody(mechanism, id))
# end
# function ∂g∂ʳposb(mechanism, eqc::EqualityConstraint{T,N,Nc,Cs}, id::Integer) where {T,N,Nc,Cs<:Tuple{<:Friction}}
#     return ∂g∂ʳpos(eqc.constraints[1], getineqconstraint(mechanism, id))
# end