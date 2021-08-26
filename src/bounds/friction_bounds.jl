mutable struct FrictionBound{T,N} <: Bound{T,N}    

    function FrictionBound(frictionid, impactid)
        new{Float64,2}(), frictionid, impactid
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

## Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳposa(ffl::FrictionBound{T,N}, λ::AbstractVector, γ::AbstractVector) where {T,N}
    return sones(T,N)
end
@inline function ∂g∂ʳposb(fvl::FrictionBound{T,N}, λ::AbstractVector, γ::AbstractVector) where {T,N}
    return szeros(T,1,1)
end
@inline function ∂g∂ʳpos(fvl::FrictionBound, λ::AbstractVector)
    return -I
end

@inline function ∂g∂ʳvela(ffl::FrictionBound{T,N}, λ::AbstractVector, γ::AbstractVector) where {T,N}
    return szeros(T,1,1)
end
@inline function ∂g∂ʳvelb(ffl::FrictionBound{T,N}, λ::AbstractVector, γ::AbstractVector) where {T,N}
    return sones(T,1,1)*ffl.cf
end

