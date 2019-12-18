mutable struct Constraint{T,N,Nc,Cs} <: Node{T,N}
    id::Int64

    constr::Cs
    pid::Int64
    linkids::SVector{Nc,Int64}

    s0::SVector{N,T}
    s1::SVector{N,T}

    function Constraint(jointdata...)
        T = jointdata[1][1].T
        linkids = Vector{Int64}(undef,0)
        constr = Vector{Joint{T}}(undef,0)
        pid = jointdata[1][2]
        N = 0
        for set in jointdata
            push!(constr,set[1])
            @assert set[2] == pid
            push!(linkids,set[3])
            N+=set[1].Nc
        end
        constr = Tuple(constr)

        Nc = length(constr)

        id = getGlobalID()

        s0 = @SVector zeros(T,N)
        s1 = @SVector zeros(T,N)

        new{T,N,Nc,typeof(constr)}(id,constr,pid,linkids,s0,s1)
    end
end

@generated function g(robot,C::Constraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(C.constr[$i],getlink(robot,C.pid),getlink(robot,C.linkids[$i]))) for i=1:Nc]
    :(vcat($(vec...)))
end

function ∂g∂pos(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    if id == C.pid
        return ∂g∂posa(robot,C,id)
    else
        return ∂g∂posb(robot,C,id)
    end
end

@generated function ∂g∂posa(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂posa(C.constr[$i],getlink(robot,id),getlink(robot,C.linkids[$i]))) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂posb(C.constr[$i],getlink(robot,C.pid),getlink(robot,id))) for i=1:Nc]
    return :(vcat($(vec...)))
end

function ∂g∂vel(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    if id == C.pid
        return ∂g∂vela(robot,C,id)
    else
        return ∂g∂velb(robot,C,id)
    end
end

@generated function ∂g∂vela(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂vela(C.constr[$i],getlink(robot,id),getlink(robot,C.linkids[$i]))) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂velb(C.constr[$i],getlink(robot,C.pid),getlink(robot,id))) for i=1:Nc]
    return :(vcat($(vec...)))
end
