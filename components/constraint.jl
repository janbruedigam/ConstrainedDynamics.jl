mutable struct Constraint{T,N,Nc,Nl,Cs} <: Node{T,N}
    id::Int64
    linkids::SVector{Nl,Int64}

    constr::Cs
    parentlink::Int64
    links::SVector{Nc,Int64}

    s0::SVector{N,T}
    s1::SVector{N,T}

    function Constraint(jointdata...)
        T = jointdata[1][1].T
        constr = Vector{Joint{T}}(undef,0)
        links = Vector{Link{T}}(undef,0)
        parentlink = jointdata[1][2]
        N = 0
        for set in jointdata
            push!(constr,set[1])
            @assert set[2] == parentlink
            push!(links,set[3])
            N+=set[1].Nc
        end
        constr = Tuple(constr)
        links = Tuple(links)

        Nc = length(constr)
        Nl = length(links)

        id = getGlobalID()
        linkids = [parentlink.id]
        for i=1:length(links)
            push!(linkids,links[i].id)
        end
        links = linkids[2:end]
        parentlink = parentlink.id
        linkids = unique(linkids)

        s0 = @SVector zeros(T,N)
        s1 = @SVector zeros(T,N)

        new{T,N,Nc,Nl,typeof(constr)}(id,linkids,constr,parentlink,links,s0,s1)
    end
end

@generated function g(robot,C::Constraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(C.constr[$i],getlink(robot,C.parentlink),getlink(robot,C.links[$i]))) for i=1:Nc]
    :(vcat($(vec...)))
end

function ∂g∂pos(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    if id == C.parentlink
        return ∂g∂posa(robot,C,id)
    else
        return ∂g∂posb(robot,C,id)
    end
end

@generated function ∂g∂posa(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂posa(C.constr[$i],getlink(robot,id),getlink(robot,C.links[$i]))) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂posb(C.constr[$i],getlink(robot,C.parentlink),getlink(robot,id))) for i=1:Nc]
    return :(vcat($(vec...)))
end

function ∂g∂vel(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    if id == C.parentlink
        return ∂g∂vela(robot,C,id)
    else
        return ∂g∂velb(robot,C,id)
    end
end

@generated function ∂g∂vela(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂vela(C.constr[$i],getlink(robot,id),getlink(robot,C.links[$i]))) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(robot,C::Constraint{T,N,Nc},id::Int64) where {T,N,Nc}
    vec = [:(∂g∂velb(C.constr[$i],getlink(robot,C.parentlink),getlink(robot,id))) for i=1:Nc]
    return :(vcat($(vec...)))
end
