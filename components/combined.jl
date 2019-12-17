#TODO vectorize constraints and links
struct Combined{T,N,Nc,Nl,Cs,L,Ls} <: Constraint{T,N,Nc,Nl}
    id::Int64
    linkids::SVector{Nl,Int64}

    constr::Cs
    parentlink::L
    links::Ls

    s0::SVector{N,T}
    s1::SVector{N,T}

    function Combined(jointdata...)
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

@generated function g(C::Combined{T,N,Nc}) where {T,N,Nc}
    vec = [:(g(C.constr[$i],C.parentlink,C.links[$i])) for i=1:Nc]
    :(vcat($(vec...)))
end

# @generated function g(C::Combined,L::Link)
#
# end

function ∂g∂pos(C::Combined,L::Link)
    id = L.id
    ids = C.linkids
    if id == ids[1]
        return [∂g∂posa(C.constr1,L,C.link2);∂g∂posa(C.constr2,L,C.link2);∂g∂posa(C.constr3,L,C.link3)]
    elseif id == ids[2]
        return [∂g∂posb(C.constr1,C.link1,L);∂g∂posb(C.constr2,C.link1,L); ∂g∂posb(C.constr3)] #[...,...,0]
    elseif id == ids[3]
        return [∂g∂posb(C.constr1);∂g∂posb(C.constr2); ∂g∂posb(C.constr3,C.link1,L)] #[0,0,...]
    else
        return [∂g∂posb(C.constr1);∂g∂posb(C.constr2);∂g∂posb(C.constr3)] #[0,0,0]
    end
end

function ∂g∂vel(C::Combined,L::Link)
    id = L.id
    ids = C.linkids
    if id == ids[1]
        return [∂g∂vela(C.constr1,L,C.link2);∂g∂vela(C.constr2,L,C.link2);∂g∂vela(C.constr3,L,C.link3)]
    elseif id == ids[2]
        return [∂g∂velb(C.constr1,C.link1,L);∂g∂velb(C.constr2,C.link1,L); ∂g∂velb(C.constr3)] #[...,...,0]
    elseif id == ids[3]
        return [∂g∂velb(C.constr1);∂g∂velb(C.constr2); ∂g∂velb(C.constr3,C.link1,L)] #[0,0,...]
    else
        return [∂g∂velb(C.constr1);∂g∂velb(C.constr2);∂g∂velb(C.constr3)] #[0,0,0]
    end
end
