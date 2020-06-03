mutable struct EqualityConstraint{T,N,Nc,Cs} <: AbstractConstraint{T,N}
    id::Int64
    name::String

    constraints::Cs
    pid::Union{Int64,Nothing}
    childids::SVector{Nc,Int64}
    inds::SVector{Nc,SVector{2,Int64}} # indices for minimal coordinates, assumes joints

    λsol::Vector{SVector{N,T}}

    function EqualityConstraint(data...; name::String="")
        jointdata = Tuple{Joint,Int64,Int64}[]
        for info in data
            if typeof(info[1]) <: Joint
                push!(jointdata, info)
            else
                for subinfo in info
                    push!(jointdata, subinfo)
                end
            end
        end

        T = getT(jointdata[1][1])# .T

        pid = jointdata[1][2]
        childids = Int64[]
        constraints = Joint{T}[]
        inds = Vector{Int64}[]
        N = 0
        for set in jointdata
            push!(constraints, set[1])
            @assert set[2] == pid
            push!(childids, set[3])

            Nset = getN(set[1])
            if isempty(inds)
                push!(inds, [1;3-Nset])
            else
                push!(inds, [last(inds)[2]+1;last(inds)[2]+3-Nset])
            end
            N += Nset
        end
        constraints = Tuple(constraints)
        Nc = length(constraints)
        

        λsol = [zeros(T, N) for i=1:2]

        new{T,N,Nc,typeof(constraints)}(getGlobalID(), name, constraints, pid, childids, inds, λsol)
    end
end

Base.length(::EqualityConstraint{T,N}) where {T,N} = N

# TODO make zero alloc
# TODO currently assumed constraints are in order and only joints which is the case unless very low level constraint setting
function setPosition!(mechanism, eqc::EqualityConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    @assert length(xθ)==3*Nc-N
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.pid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.childids[i])
        Δx = getPositionDelta(eqc.constraints[i], body1, body2, xθ[SUnitRange(eqc.inds[i][1],eqc.inds[i][2])]) 
        Δq = getPositionDelta(eqc.constraints[i+1], body1, body2, xθ[SUnitRange(eqc.inds[i+1][1],eqc.inds[i+1][2])])
        
        p1, p2 = eqc.constraints[i].vertices
        setPosition!(body1, body2; p1 = p1, p2 = p2, Δx = Δx, Δq = Δq)
    end
end

# TODO make zero alloc
# TODO currently assumed constraints are in order and only joints which is the case unless very low level constraint setting
function setVelocity!(mechanism, eqc::EqualityConstraint{T,N,Nc}, vω) where {T,N,Nc}
    @assert length(vω)==3*Nc-N
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.pid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.childids[i])
        Δv = getVelocityDelta(eqc.constraints[i], body1, body2, vω[SUnitRange(eqc.inds[i][1],eqc.inds[i][2])])
        Δω = getVelocityDelta(eqc.constraints[i+1], body1, body2, vω[SUnitRange(eqc.inds[i+1][1],eqc.inds[i+1][2])])
        
        p1, p2 = eqc.constraints[i].vertices
        setVelocity!(body1, body2; p1 = p1, p2 = p2, Δv = Δv, Δω = Δω)
    end
end

# TODO make zero alloc
function setForce!(mechanism, eqc::EqualityConstraint{T,N,Nc}, Fτ::AbstractVector{T}) where {T,N,Nc}
    @assert length(Fτ)==3*Nc-N
    for i = 1:Nc
        setForce!(eqc.constraints[i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.childids[i]),  Fτ[SUnitRange(eqc.inds[i][1],eqc.inds[i][2])])
    end
    return
end

@generated function minimalCoordinates(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimalCoordinates(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    :(svcat($(vec...)))
end

@inline function GtλTof!(mechanism, body::Body, eqc::EqualityConstraint)
    body.f -= zerodimstaticadjoint(∂g∂pos(mechanism, eqc, body.id)) * eqc.λsol[2]
    return
end

@generated function g(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.childids[$i]), mechanism.Δt)) for i = 1:Nc]
    :(svcat($(vec...)))
end

@inline function ∂g∂pos(mechanism, eqc::EqualityConstraint, id::Integer)
    id == eqc.pid ? ∂g∂posa(mechanism, eqc, id) : ∂g∂posb(mechanism, eqc, id)
end

@inline function ∂g∂vel(mechanism, eqc::EqualityConstraint, id::Integer)
    id == eqc.pid ? ∂g∂vela(mechanism, eqc, id) : ∂g∂velb(mechanism, eqc, id)
end

@inline function ∂g∂con(mechanism, eqc::EqualityConstraint, id::Integer)
    ∂g∂con(mechanism, eqc)
end

@generated function ∂g∂posa(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂posa(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.childids[$i]))) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂posb(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, id))) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂vela(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂vela(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.childids[$i]), mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Integer) where {T,N,Nc}
    vec = [:(∂g∂velb(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, id), mechanism.Δt)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂con(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(∂g∂con(eqc.constraints[$i])) for i = 1:Nc]
    return :(vcat($(vec...)))
end