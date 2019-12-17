mutable struct Constraint{T,N,Nc,Nl,Cs,L,Ls} <: Node{T,N}
    id::Int64
    linkids::SVector{Nl,Int64}

    constr::Cs
    parentlink::L
    links::Ls

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
        linkids = unique(linkids)

        s0 = @SVector zeros(T,N)
        s1 = @SVector zeros(T,N)

        new{T,N,Nc,Nl,typeof(constr),typeof(parentlink),typeof(links)}(id,linkids,constr,parentlink,links,s0,s1)
    end
end

@generated function g(C::Constraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(C.constr[$i],C.parentlink,C.links[$i])) for i=1:Nc]
    :(vcat($(vec...)))
end

function ∂g∂pos(C::Constraint{T,N,Nc},L::Link) where {T,N,Nc}
    id = L.id
    if id == C.parentlink.id
        return ∂g∂posa(C,L)
    else
        return ∂g∂posb(C,L)
    end
end

@generated function ∂g∂posa(C::Constraint{T,N,Nc},L::Link) where {T,N,Nc}
    vec = [:(∂g∂posa(C.constr[$i],L,C.links[$i])) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂posb(C::Constraint{T,N,Nc},L::Link) where {T,N,Nc}
    vec = [:(∂g∂posb(C.constr[$i],C.parentlink,L)) for i=1:Nc]
    return :(vcat($(vec...)))
end

function ∂g∂vel(C::Constraint{T,N,Nc},L::Link) where {T,N,Nc}
    id = L.id
    if id == C.parentlink.id
        return ∂g∂vela(C,L)
    else
        return ∂g∂velb(C,L)
    end
end

@generated function ∂g∂vela(C::Constraint{T,N,Nc},L::Link) where {T,N,Nc}
    vec = [:(∂g∂vela(C.constr[$i],L,C.links[$i])) for i=1:Nc]
    return :(vcat($(vec...)))
end

@generated function ∂g∂velb(C::Constraint{T,N,Nc},L::Link) where {T,N,Nc}
    vec = [:(∂g∂velb(C.constr[$i],C.parentlink,L)) for i=1:Nc]
    return :(vcat($(vec...)))
end
