mutable struct EqualityConstraint{T,N,Nc,Cs} <: AbstractConstraint{T,N}
    id::Int64

    constraints::Cs
    pid::Union{Int64,Nothing}
    bodyids::SVector{Nc,Int64}

    s0::SVector{N,T}
    s1::SVector{N,T}

    function EqualityConstraint(data...)
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

        new{T,N,Nc,typeof(constraints)}(getGlobalID(), constraints, pid, bodyids, s0, s1)
    end
end

Base.length(::EqualityConstraint{T,N}) where {T,N} = N

function setForce!(Fτ, eqc::EqualityConstraint{T,N,Nc}, mechanism; K=mechanism.No) where {T,N,Nc}
    for i = 1:Nc
        setForce!(eqc.constraints[i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.bodyids[i]), Fτ[i], K)
    end
end

@generated function minimalCoordinates(eqc::EqualityConstraint{T,N,Nc}, mechanism; K=mechanism.No) where {T,N,Nc}
    vec = [:(minimalCoordinates(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.bodyids[$i]), K)) for i = 1:Nc]
    :(svcat($(vec...)))
end

@generated function g(eqc::EqualityConstraint{T,N,Nc}, mechanism) where {T,N,Nc}
    vec = [:(g(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, eqc.bodyids[$i]), mechanism.Δt, mechanism.No)) for i = 1:Nc]
    :(svcat($(vec...)))
end

@inline function ∂g∂pos(eqc::EqualityConstraint, id::Int64, mechanism)
    id == eqc.pid ? ∂g∂posa(eqc, id, mechanism) : ∂g∂posb(eqc, id, mechanism)
end

@inline function ∂g∂vel(eqc::EqualityConstraint, id::Int64, mechanism)
    id == eqc.pid ? ∂g∂vela(eqc, id, mechanism) : ∂g∂velb(eqc, id, mechanism)
end

@generated function ∂g∂posa(eqc::EqualityConstraint{T,N,Nc}, id::Int64, mechanism) where {T,N,Nc}
    vec = [:(∂g∂posa(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.bodyids[$i]), mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(eqc::EqualityConstraint{T,N,Nc}, id::Int64, mechanism) where {T,N,Nc}
    vec = [:(∂g∂posb(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, id), mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂vela(eqc::EqualityConstraint{T,N,Nc}, id::Int64, mechanism) where {T,N,Nc}
    vec = [:(∂g∂vela(eqc.constraints[$i], getbody(mechanism, id), getbody(mechanism, eqc.bodyids[$i]), mechanism.Δt, mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(eqc::EqualityConstraint{T,N,Nc}, id::Int64, mechanism) where {T,N,Nc}
    vec = [:(∂g∂velb(eqc.constraints[$i], getbody(mechanism, eqc.pid), getbody(mechanism, id), mechanism.Δt, mechanism.No)) for i = 1:Nc]
    return :(vcat($(vec...)))
end
