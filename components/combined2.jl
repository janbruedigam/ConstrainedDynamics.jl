#TODO vectorize constraints and links
struct Combined2{T,Nc,Nl,C1,C2,L1,L2} <: Constraint{T,Nc,Nl}
    id::Int64
    linkids::SVector{Nl,Int64}

    constr1::C1
    constr2::C2
    link1::L1
    link2::L2

    s0::SVector{Nc,T}
    s1::SVector{Nc,T}

    data::NodeData{T,Nc}

    function Combined2(c1, c2)
        constr1,l1,l2 = c1
        constr2,l3,l4 = c2
        constr = [constr1;constr2]
        links = unique([l1;l2;l3;l4])

        T = constr1.T
        Nl = length(links)
        Nc = constr1.Nc + constr2.Nc

        id = getGlobalID()
        linkids = unique([links[i].id for i=1:Nl])

        s0 = @SVector zeros(T,Nc)
        s1 = @SVector zeros(T,Nc)

        data = NodeData{T,Nc}()

        new{T,Nc,Nl,typeof(constr[1]),typeof(constr[2]),typeof(links[1]),typeof(links[2])}(id,linkids,constr...,links...,s0,s1,data)
    end
end

g(C::Combined2) = [g(C.constr1,C.link1,C.link2);g(C.constr2,C.link1,C.link2)]

function ∂g∂pos(C::Combined2,L::Link)
    id = L.id
    ids = C.linkids
    if id == ids[1]
        return [∂g∂posa(C.constr1,L,C.link2);∂g∂posa(C.constr2,L,C.link2)]
    elseif id == ids[2]
        return [∂g∂posb(C.constr1,C.link1,L);∂g∂posb(C.constr2,C.link1,L)]
    else
        return [∂g∂posb(C.constr1);∂g∂posb(C.constr2)] #[0,0]
    end
end

function ∂g∂vel(C::Combined2,L::Link)
    id = L.id
    ids = C.linkids
    if id == ids[1]
        return [∂g∂vela(C.constr1,L,C.link2);∂g∂vela(C.constr2,L,C.link2)]
    elseif id == ids[2]
        return [∂g∂velb(C.constr1,C.link1,L);∂g∂velb(C.constr2,C.link1,L)]
    else
        return [∂g∂velb(C.constr1);∂g∂velb(C.constr2)] #[0,0]
    end
end

getNc(C::Combined2) = C.constr1.Nc+C.constr2.Nc
# @generated function linkids(C::Combined2{T,Nc,Nc²,Nl}) where {T,Nc,Nc²,Nl}
#     ids = [:(C.links[$i].data.id) for i=1:Nl]
#     :(SVector{Nl,Int64}($(ids...)))
# end
