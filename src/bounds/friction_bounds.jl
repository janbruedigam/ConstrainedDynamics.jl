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
@inline g(ffl::FrictionBound, β::AbstractVector, γ::AbstractVector) = ffl.cf*γ - SVector{1,Float64}(sum(β))
@inline g(fvl::BetaBound, β::AbstractVector) = β

## Derivatives 
@inline function ∂g∂beta(::Friction{T,N}, ::FrictionBound) where {T,N}
    return [szeros(T,N) sones(T,N)]
end
@inline function ∂g∂beta(::Friction{T,N}, ::BetaBound) where {T,N}
    return [szeros(T,N,N) -I]
end

@inline function ∂gab∂ʳba(ffl::FrictionBound{T}, impact::Impact) where T
    SA{T}[0 0;0 ffl.cf], szeros(T,2,2)
end
