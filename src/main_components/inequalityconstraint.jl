mutable struct InequalityConstraint{T,N,C} <: AbstractConstraint{T,N}
    id::Int64
    name::String

    constraint::C
    parentid::Int64
    # childid::Int64

    ssol::Vector{SVector{N,T}}
    γsol::Vector{SVector{N,T}}
    

    function InequalityConstraint(data; name::String="")
        bound, id = data
        T = getT(bound)

        parentid = id
        # childids = Int64[]
        constraint = bound
        N = length(constraint)

        ssol = [ones(T, N) for i=1:2]
        γsol = [ones(T, N) for i=1:2]

        new{T,N,typeof(constraint)}(getGlobalID(), name, constraint, parentid, ssol, γsol)
    end
end


function resetVars!(ineqc::InequalityConstraint{T,N}) where {T,N}
    ineqc.ssol[1] = sones(T, N)
    ineqc.ssol[2] = sones(T, N)
    ineqc.γsol[1] = sones(T, N)
    ineqc.γsol[2] = sones(T, N)

    return 
end

@inline function NtγTof!(mechanism, body::Body, ineqc::InequalityConstraint{T,N}) where {T,N}
    body.state.d -= ∂g∂ʳpos(mechanism, ineqc, body.id)' * ineqc.γsol[2]
    return
end

function g(mechanism, ineqc::InequalityConstraint)
    return g(ineqc.constraint, getbody(mechanism, ineqc.parentid), mechanism.Δt)
end

function gs(mechanism, ineqc::InequalityConstraint)
    return g(ineqc.constraint, getbody(mechanism, ineqc.parentid), mechanism.Δt) - ineqc.ssol[2]
end

function complementarity(mechanism, ineqc::InequalityConstraint)
    return ineqc.γsol[2] .* ineqc.ssol[2]
end

function complementarityμ(mechanism, ineqc::InequalityConstraint)
    return ineqc.γsol[2] .* ineqc.ssol[2] .- mechanism.μ
end

function ∂g∂ʳpos(mechanism, ineqc::InequalityConstraint, id::Integer)
    return ∂g∂ʳpos(ineqc.constraint, getbody(mechanism, id))
end

function ∂g∂ʳvel(mechanism, ineqc::InequalityConstraint, id::Integer)
    return ∂g∂ʳvel(ineqc.constraint, getbody(mechanism, id), mechanism.Δt)
end