# TODO currently just single parent and child
mutable struct InequalityConstraint{T,N,Nc,Cs,N½} <: AbstractConstraint{T,N}
    id::Int64
    name::String

    # Currently only single constraint and child
    constraints::Cs
    parentid::Int64
    childids::SVector{1,Union{Int64,Nothing}}

    ssol::Vector{SVector{N½,T}}
    γsol::Vector{SVector{N½,T}}
    

    function InequalityConstraint(data; name::String="")
        bound, parentid, childid = data
        T = getT(bound)

        childids = [childid]
        constraint = Tuple([bound])
        N = length(constraint[1])
        N½ = Int64(N/2)

        ssol = [ones(T, N½) for i=1:2]
        γsol = [ones(T, N½) for i=1:2]

        new{T,N,1,typeof(constraint),N½}(getGlobalID(), name, constraint, parentid, childids, ssol, γsol)
    end
end


function resetVars!(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    ineqc.ssol[1] = sones(T, N½)
    ineqc.ssol[2] = sones(T, N½)
    ineqc.γsol[1] = sones(T, N½)
    ineqc.γsol[2] = sones(T, N½)

    return 
end

@inline function constraintForceMapping!(mechanism, body::Body, ineqc::InequalityConstraint)
    body.state.d -= ∂g∂ʳpos(mechanism, ineqc, body)' * ineqc.γsol[2]
    return
end

# function g(mechanism, ineqc::InequalityConstraint)
#     childid = ineqc.childids[1]
#     if childid === nothing
#         return g(ineqc.constraints[1], getcomponent(mechanism, ineqc.parentid), mechanism.Δt)
#     else
#         return g(ineqc.constraints[1], getfriction(mechanism, ineqc.parentid), getineqconstraint(mechanism, childid))
#     end
# end

function gs(mechanism, ineqc::InequalityConstraint)
    return g(mechanism, ineqc) - ineqc.ssol[2]
end

function complementarity(mechanism, ineqc::InequalityConstraint)
    return ineqc.γsol[2] .* ineqc.ssol[2]
end
function complementarityμ(mechanism, ineqc::InequalityConstraint)
    return ineqc.γsol[2] .* ineqc.ssol[2] .- mechanism.μ
end

@inline function ∂gab∂ʳba(mechanism, body::Body, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z;-∂g∂ʳpos(mechanism, ineqc, body)]', [Z;∂g∂ʳvel(mechanism, ineqc, body)]
end
@inline function ∂gab∂ʳba(mechanism, ineqc1::InequalityConstraint, ineqc2::InequalityConstraint)
    G1, G2 = ∂gab∂ʳba(ineqc1.constraints[1], ineqc2.constraints[1])

    return G1, G2
end

function ∂g∂ʳposa(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳposa(ineqc.constraints[1], body, nothing)
end
# function ∂g∂ʳposb(mechanism, ineqc::InequalityConstraint, body::Body)
#     return ∂g∂ʳposb(ineqc.constraints[1], body, nothing)
# end

function ∂g∂ʳvela(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳvela(ineqc.constraints[1], body, nothing, mechanism.Δt)
end
# function ∂g∂ʳvelb(mechanism, ineqc::InequalityConstraint, body::Body)
#     return ∂g∂ʳvelb(ineqc.constraints[1], body, nothing, mechanism.Δt)
# end
