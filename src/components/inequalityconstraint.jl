mutable struct InequalityConstraint{T,N,Cs} <: AbstractConstraint{T,N}
    id::Int64
    name::String
    active::Bool

    constraints::Cs
    parentid::Int64
    # childid::Int64

    ssol::Vector{SVector{N,T}}
    γsol::Vector{SVector{N,T}}
    

    function InequalityConstraint(data...; name::String="")
        bounddata = Tuple{Bound,Int64}[]
        for info in data
            if typeof(info[1]) <: Bound
                push!(bounddata, info)
            else
                for subinfo in info
                    push!(bounddata, subinfo)
                end
            end
        end

        T = getT(bounddata[1][1])

        parentid = bounddata[1][2]
        childids = Int64[]
        constraints = Bound{T}[]
        N = 0
        for set in bounddata
            push!(constraints, set[1])
            @assert set[2] == parentid
            N += 1 # getN(set[1])
        end
        constraints = Tuple(constraints)
        # Nc = length(constraints)

        ssol = [ones(T, N) for i=1:2]
        γsol = [ones(T, N) for i=1:2]

        new{T,N,typeof(constraints)}(getGlobalID(), name, true, constraints, parentid, ssol, γsol)
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
    if isactive(ineqc)
        state = body.state
        state.d -= ∂g∂pos(mechanism, ineqc, body)' * ineqc.γsol[2]
        for i=1:N
            state.d -= additionalforce(ineqc.constraints[i])
        end
    end
    return
end

function g(mechanism, ineqc::InequalityConstraint{T,1}) where {T}
    return g(ineqc.constraints[1], getbody(mechanism, ineqc.parentid), mechanism.Δt)
end

@generated function g(mechanism, ineqc::InequalityConstraint{T,N}) where {T,N}
    vec = [:(g(ineqc.constraints[$i], getbody(mechanism, ineqc.parentid), mechanism.Δt)) for i = 1:N]
    return :(SVector{N,T}($(vec...)))
end

function gs(mechanism, ineqc::InequalityConstraint{T,1}) where {T}
    return g(ineqc.constraints[1], getbody(mechanism, ineqc.parentid), mechanism.Δt) - ineqc.ssol[2][1]
end

@generated function gs(mechanism, ineqc::InequalityConstraint{T,N}) where {T,N}
    vec = [:(g(ineqc.constraints[$i], getbody(mechanism, ineqc.parentid), mechanism.Δt) - ineqc.ssol[2][$i]) for i = 1:N]
    return :(SVector{N,T}($(vec...)))
end

function h(ineqc::InequalityConstraint)
    return ineqc.ssol[2] .* ineqc.γsol[2]
end

function hμ(ineqc::InequalityConstraint{T}, μ) where T
    return ineqc.ssol[2] .* ineqc.γsol[2] .- μ
end


function schurf(mechanism, ineqc::InequalityConstraint{T,N}, body) where {T,N}
    val = szeros(T, 6)
    for i = 1:N
        val += schurf(ineqc, ineqc.constraints[i], i, body, mechanism.μ, mechanism.Δt)
    end
    return val
end

function schurD(ineqc::InequalityConstraint{T,N}, body, Δt) where {T,N}
    val = szeros(T, 6, 6)
    for i = 1:N
        val += schurD(ineqc, ineqc.constraints[i], i, body, Δt)
    end
    return val
end

@generated function ∂g∂pos(mechanism, ineqc::InequalityConstraint{T,N}, body) where {T,N}
    vec = [:(∂g∂pos(ineqc.constraints[$i])) for i = 1:N]
    return :(vcat($(vec...)))
end

@generated function ∂g∂vel(mechanism, ineqc::InequalityConstraint{T,N}, body) where {T,N}
    vec = [:(∂g∂vel(ineqc.constraints[$i], mechanism.Δt)) for i = 1:N]
    return :(vcat($(vec...)))
end

function calcFrictionForce!(mechanism, ineqc::InequalityConstraint{T,N}) where {T,N}
    for i = 1:N
        constraint = ineqc.constraints[i]
        if typeof(constraint) <: Friction
            calcFrictionForce!(mechanism, ineqc, constraint, i, getbody(mechanism, ineqc.parentid))
        end
    end
    return
end