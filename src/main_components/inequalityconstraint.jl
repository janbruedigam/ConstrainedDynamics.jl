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

        # childids = Int64[]
        constraint = Tuple([bound])
        N = length(constraint[1])
        N½ = Int64(N/2)

        ssol = [ones(T, N½) for i=1:2]
        γsol = [ones(T, N½) for i=1:2]

        new{T,N,1,typeof(constraint),N½}(getGlobalID(), name, constraint, parentid, [childid], ssol, γsol)
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

function g(mechanism, ineqc::InequalityConstraint)
    return g(ineqc.constraints[1], getbody(mechanism, ineqc.parentid), mechanism.Δt)
end

#TODO how to do allocation free dispatch better?
# function g(mechanism::AbstractMechanism{T,Nn,Nb,Ne}, ineqc::InequalityConstraint) where {T,Nn,Nb,Ne}
#     parentid = ineqc.parentid
#     childid = ineqc.childids[1]
#     if parentid <= Ne
#         if childid === nothing
#             return g(ineqc.constraints[1], geteqconstraint(mechanism, parentid), mechanism.Δt)
#         elseif childid <= Ne
#             # return g(ineqc.constraints[1], geteqconstraint(mechanism, parentid), geteqconstraint(mechanism, childid), mechanism.Δt)
#         elseif childid <= Ne+Nb
#             # return g(ineqc.constraints[1], geteqconstraint(mechanism, parentid), getbody(mechanism, childid), mechanism.Δt)
#         else
#             return g(ineqc.constraints[1], geteqconstraint(mechanism, parentid), getineqconstraint(mechanism, childid), mechanism.Δt)
#         end
#     elseif parentid <= Ne+Nb
#         if childid === nothing
#             return g(ineqc.constraints[1], getbody(mechanism, parentid), mechanism.Δt)
#         elseif childid <= Ne
#             # return g(ineqc.constraints[1], getbody(mechanism, parentid), geteqconstraint(mechanism, childid), mechanism.Δt)
#         elseif childid <= Ne+Nb
#             # return g(ineqc.constraints[1], getbody(mechanism, parentid), getbody(mechanism, childid), mechanism.Δt)
#         else
#             # return g(ineqc.constraints[1], getbody(mechanism, parentid), getineqconstraint(mechanism, childid), mechanism.Δt)
#         end
#     # else
#     #     if childid === nothing
#     #         # return g(ineqc.constraints[1], getineqconstraint(mechanism, parentid), mechanism.Δt)
#     #     elseif childid <= Ne
#     #         # return g(ineqc.constraints[1], getineqconstraint(mechanism, parentid), geteqconstraint(mechanism, childid), mechanism.Δt)
#     #     elseif childid <= Ne+Nb
#     #         # return g(ineqc.constraints[1], getineqconstraint(mechanism, parentid), getbody(mechanism, childid), mechanism.Δt)
#     #     else
#     #         # return g(ineqc.constraints[1], getineqconstraint(mechanism, parentid), getineqconstraint(mechanism, childid), mechanism.Δt)
#     #     end
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

# # TODO these are currently specialized to contact and friction constraints
# function ∂g∂ʳposa(mechanism::AbstractMechanism{T,Nn,Nb,Ne}, ineqc::InequalityConstraint, id::Integer) where {T,Nn,Nb,Ne}
#     childid = ineqc.childids[1]
#     if id <= Ne
#         if childid === nothing
#             return ∂g∂ʳposa(ineqc.constraints[1], geteqconstraint(mechanism, id), childid)
#         elseif childid <= Ne
#             # return ∂g∂ʳposa(ineqc.constraints[1], geteqconstraint(mechanism, id), geteqconstraint(mechanism, childid), childid)
#         elseif childid <= Ne+Nb
#             # return ∂g∂ʳposa(ineqc.constraints[1], geteqconstraint(mechanism, id), getbody(mechanism, childid), childid)
#         else
#             return ∂g∂ʳposa(ineqc.constraints[1], geteqconstraint(mechanism, id), getineqconstraint(mechanism, childid), childid)
#         end
#     elseif id <= Ne+Nb
#         if childid === nothing
#             return ∂g∂ʳposa(ineqc.constraints[1], getbody(mechanism, id), childid)
#         elseif childid <= Ne
#             # return ∂g∂ʳposa(ineqc.constraints[1], getbody(mechanism, id), geteqconstraint(mechanism, childid), childid)
#         elseif childid <= Ne+Nb
#             # return ∂g∂ʳposa(ineqc.constraints[1], getbody(mechanism, id), getbody(mechanism, childid), childid)
#         else
#             # return ∂g∂ʳposa(ineqc.constraints[1], getbody(mechanism, id), getineqconstraint(mechanism, childid), childid)
#         end
#     # else
#     #     if childid === nothing
#     #         # return ∂g∂ʳposa(ineqc.constraints[1], getineqconstraint(mechanism, id), childid)
#     #     elseif childid <= Ne
#     #         # return ∂g∂ʳposa(ineqc.constraints[1], getineqconstraint(mechanism, id), geteqconstraint(mechanism, childid), childid)
#     #     elseif childid <= Ne+Nb
#     #         # return ∂g∂ʳposa(ineqc.constraints[1], getineqconstraint(mechanism, id), getbody(mechanism, childid), childid)
#     #     else
#     #         # return ∂g∂ʳposa(ineqc.constraints[1], getineqconstraint(mechanism, id), getineqconstraint(mechanism, childid), childid)
#     #     end
#     end
# end
# function ∂g∂ʳposb(mechanism::AbstractMechanism{T,Nn,Nb,Ne}, ineqc::InequalityConstraint, id::Integer) where {T,Nn,Nb,Ne}
#     parentid = ineqc.parentid
#     if parentid <= Ne
#         if id === nothing
#             # return ∂g∂ʳposb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), id)
#         elseif id <= Ne
#             # return ∂g∂ʳposb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), geteqconstraint(mechanism, id), id)
#         elseif id <= Ne+Nb
#             # return ∂g∂ʳposb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), getbody(mechanism, id), id)
#         else
#             return ∂g∂ʳposb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), getineqconstraint(mechanism, id), id)
#         end
#     # elseif parentid <= Ne+Nb
#     #     if id === nothing
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getbody(mechanism, parentid), id)
#     #     elseif id <= Ne
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getbody(mechanism, parentid), geteqconstraint(mechanism, id), id)
#     #     elseif id <= Ne+Nb
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getbody(mechanism, parentid), getbody(mechanism, id), id)
#     #     else
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getbody(mechanism, parentid), getineqconstraint(mechanism, id), id)
#     #     end
#     # else
#     #     if id === nothing
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), id)
#     #     elseif id <= Ne
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), geteqconstraint(mechanism, id), id)
#     #     elseif id <= Ne+Nb
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), getbody(mechanism, id), id)
#     #     else
#     #         # return ∂g∂ʳposb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), getineqconstraint(mechanism, id), id)
#     #     end
#     end
# end

# function ∂g∂ʳvela(mechanism::AbstractMechanism{T,Nn,Nb,Ne}, ineqc::InequalityConstraint, id::Integer) where {T,Nn,Nb,Ne}
#     childid = ineqc.childids[1]
#     if id <= Ne
#         if childid === nothing
#             return ∂g∂ʳvela(ineqc.constraints[1], geteqconstraint(mechanism, id), childid, mechanism.Δt)
#         elseif childid <= Ne
#             # return ∂g∂ʳvela(ineqc.constraints[1], geteqconstraint(mechanism, id), geteqconstraint(mechanism, childid), childid, mechanism.Δt)
#         elseif childid <= Ne+Nb
#             # return ∂g∂ʳvela(ineqc.constraints[1], geteqconstraint(mechanism, id), getbody(mechanism, childid), childid, mechanism.Δt)
#         else
#             return ∂g∂ʳvela(ineqc.constraints[1], geteqconstraint(mechanism, id), getineqconstraint(mechanism, childid), childid, mechanism.Δt)
#         end
#     elseif id <= Ne+Nb
#         if childid === nothing
#             return ∂g∂ʳvela(ineqc.constraints[1], getbody(mechanism, id), childid, mechanism.Δt)
#         elseif childid <= Ne
#             # return ∂g∂ʳvela(ineqc.constraints[1], getbody(mechanism, id), geteqconstraint(mechanism, childid), childid, mechanism.Δt)
#         elseif childid <= Ne+Nb
#             # return ∂g∂ʳvela(ineqc.constraints[1], getbody(mechanism, id), getbody(mechanism, childid), childid, mechanism.Δt)
#         else
#             # return ∂g∂ʳvela(ineqc.constraints[1], getbody(mechanism, id), getineqconstraint(mechanism, childid), childid, mechanism.Δt)
#         end
#     # else
#     #     if childid === nothing
#     #         # return ∂g∂ʳvela(ineqc.constraints[1], getineqconstraint(mechanism, id), childid, mechanism.Δt)
#     #     elseif childid <= Ne
#     #         # return ∂g∂ʳvela(ineqc.constraints[1], getineqconstraint(mechanism, id), geteqconstraint(mechanism, childid), childid, mechanism.Δt)
#     #     elseif childid <= Ne+Nb
#     #         # return ∂g∂ʳvela(ineqc.constraints[1], getineqconstraint(mechanism, id), getbody(mechanism, childid), childid, mechanism.Δt)
#     #     else
#     #         # return ∂g∂ʳvela(ineqc.constraints[1], getineqconstraint(mechanism, id), getineqconstraint(mechanism, childid), childid, mechanism.Δt)
#     #     end
#     end
# end
# function ∂g∂ʳvelb(mechanism::AbstractMechanism{T,Nn,Nb,Ne}, ineqc::InequalityConstraint, id::Integer) where {T,Nn,Nb,Ne}
#     parentid = ineqc.parentid
#     if parentid <= Ne
#         if id === nothing
#             # return ∂g∂ʳvelb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), id, mechanism.Δt)
#         elseif id <= Ne
#             # return ∂g∂ʳvelb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), geteqconstraint(mechanism, id), id, mechanism.Δt)
#         elseif id <= Ne+Nb
#             # return ∂g∂ʳvelb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), getbody(mechanism, id), id, mechanism.Δt)
#         else
#             return ∂g∂ʳvelb(ineqc.constraints[1], geteqconstraint(mechanism, parentid), getineqconstraint(mechanism, id), id, mechanism.Δt)
#         end
#     # elseif parentid <= Ne+Nb
#     #     if id === nothing
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getbody(mechanism, parentid), id, mechanism.Δt)
#     #     elseif id <= Ne
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getbody(mechanism, parentid), geteqconstraint(mechanism, id), id, mechanism.Δt)
#     #     elseif id <= Ne+Nb
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getbody(mechanism, parentid), getbody(mechanism, id), id, mechanism.Δt)
#     #     else
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getbody(mechanism, parentid), getineqconstraint(mechanism, id), id, mechanism.Δt)
#     #     end
#     # else
#     #     if id === nothing
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), id, mechanism.Δt)
#     #     elseif id <= Ne
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), geteqconstraint(mechanism, id), id, mechanism.Δt)
#     #     elseif id <= Ne+Nb
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), getbody(mechanism, id), id, mechanism.Δt)
#     #     else
#     #         # return ∂g∂ʳvelb(ineqc.constraints[1], getineqconstraint(mechanism, parentid), getineqconstraint(mechanism, id), id, mechanism.Δt)
#     #     end
#     end
# end
