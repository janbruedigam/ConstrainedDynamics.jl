mutable struct FrictionBound{T,N} <: Bound{T,N}
    cf::T

    function FrictionBound(cf, frictionid, impactid)
        new{Float64,2}(cf), frictionid, impactid
    end
end

mutable struct BetaBound{T,N} <: Bound{T,N}

    function BetaBound(frictionid)
        new{Float64,8}(), frictionid, nothing
    end
end


### Constraints and derivatives
## Position level constraints (for dynamics)
function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{Impact{T,N}}}
    g(ineqc.constraints[1], getbody(mechanism, ineqc.parentid), mechanism.Δt)
end

g(bound::FrictionBound, fric::Friction, Δt) = g(bound, fric.βsol[2], fric.γsolref[2])
g(bound::BetaBound, fric::Friction, Δt) = g(bound, fric.βsol[2])

@inline g(ffl::FrictionBound, β::AbstractVector, γ::AbstractVector) = ffl.cf*γ - SVector{1,Float64}(sum(β))
@inline g(fvl::BetaBound, β::AbstractVector) = β

@inline function constraintForceMapping!(mechanism, fric::Friction, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{FrictionBound{T,N}}}
    fric.d = fric.d .+ ineqc.γsol[2]
    return
end
@inline function constraintForceMapping!(mechanism, fric::Friction, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{BetaBound{T,N}}}
    fric.d -= ineqc.γsol[2]
    return
end

## Derivatives 
@inline function ∂g∂beta(::Friction{T}, ::FrictionBound) where {T}
    return [szeros(T,4) sones(T,4)]
end
@inline function ∂g∂beta(::Friction{T}, ::BetaBound) where {T}
    return [szeros(T,4,4) -I]
end

@inline function ∂gab∂ʳba(ffl::FrictionBound{T}, impact::Impact) where T
    SA{T}[0 0;0 ffl.cf], szeros(T,2,2)
end