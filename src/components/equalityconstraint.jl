mutable struct EqualityConstraint{T,N,Nc,Cs} <: AbstractConstraint{T,N}
    id::Int64
    name::String

    constraints::Cs
    pid::Union{Int64,Nothing}
    bodyids::SVector{Nc,Int64}

    s0::SVector{N,T}
    s1::SVector{N,T}

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
        bodyids = Int64[]
        constraints = Joint{T}[]
        N = 0
        for set in jointdata
            push!(constraints, set[1])
            @assert set[2] == pid
            push!(bodyids, set[3])
            N += getNc(set[1])
        end
        constraints = Tuple(constraints)
        Nc = length(constraints)

        s0 = @SVector zeros(T, N)
        s1 = @SVector zeros(T, N)

        new{T,N,Nc,typeof(constraints)}(getGlobalID(), name, constraints, pid, bodyids, s0, s1)
    end
end

Base.length(::EqualityConstraint{T,N}) where {T,N} = N

function setForce!(mechanism, eqc::EqualityConstraint{T,N,Nc}, Fτ; K = mechanism.No) where {T,N,Nc}
    for i = 1:Nc
        setForce!(eqc.constraints[i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.bodyids[i]), Fτ[i], K)
    end
end

function setVelocity!(mechanism, eqc::EqualityConstraint{T,N,Nc}, vω) where {T,N,Nc}
    # TODO currently assumed constraints are in order which is the case unless very low level constraint setting and only joints
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.pid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.bodyids[i])
        Δv = getVelocityDelta(eqc.constraints[i], body1, body2, vω[i])
        Δω = getVelocityDelta(eqc.constraints[i+1], body1, body2, vω[i+1])
        
        p1, p2 = eqc.constraints[i].vertices
        setVelocity!(mechanism, body1, body2; p1 = p1, p2 = p2, Δv = Δv, Δω = Δω)
    end
end

function setPosition!(mechanism, eqc::EqualityConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    # TODO currently assumed constraints are in order which is the case unless very low level constraint setting and only joints
    n = Int64(Nc/2)
    body1 = getbody(mechanism, eqc.pid)
    for i = 1:n
        body2 = getbody(mechanism, eqc.bodyids[i])
        Δx = getPositionDelta(eqc.constraints[i], body1, body2, xθ[i])
        Δq = getPositionDelta(eqc.constraints[i+1], body1, body2, xθ[i+1])
        
        p1, p2 = eqc.constraints[i].vertices
        setPosition!(mechanism, body1, body2; p1 = p1, p2 = p2, Δx = Δx, Δq = Δq)
    end
end

@generated function minimalCoordinates(mechanism, eqc::EqualityConstraint{T,N,Nc}; K = mechanism.No) where {T,N,Nc}
    vec = [:(minimalCoordinates(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.bodyids[$i]), K)) for i = 1:Nc]
    :(svcat($(vec...)))
end

@generated function g(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.bodyids[$i]), mechanism.Δt, mechanism.No)) for i = 1:Nc]
    :(svcat($(vec...)))
end

@inline function ∂g∂pos(mechanism, eqc::EqualityConstraint, id::Int64)
    id == eqc.pid ? ∂g∂posa(mechanism, eqc, id) : ∂g∂posb(mechanism, eqc, id)
end

@inline function ∂g∂vel(mechanism, eqc::EqualityConstraint, id::Int64)
    id == eqc.pid ? ∂g∂vela(mechanism, eqc, id) : ∂g∂velb(mechanism, eqc, id)
end

@inline function ∂g∂con(mechanism, eqc::EqualityConstraint, id::Int64)
    ∂g∂con(mechanism, eqc)
end

@generated function ∂g∂posa(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Int64) where {T,N,Nc}
    vec = [:(∂g∂posa(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.bodyids[$i]), mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Int64) where {T,N,Nc}
    vec = [:(∂g∂posb(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, id), mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂vela(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Int64) where {T,N,Nc}
    vec = [:(∂g∂vela(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.bodyids[$i]), mechanism.Δt, mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(mechanism, eqc::EqualityConstraint{T,N,Nc}, id::Int64) where {T,N,Nc}
    vec = [:(∂g∂velb(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, id), mechanism.Δt, mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂con(mechanism, eqc::EqualityConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(∂g∂con(eqc.constraints[$i])) for i = 1:Nc]
    return :(vcat($(vec...)))
end