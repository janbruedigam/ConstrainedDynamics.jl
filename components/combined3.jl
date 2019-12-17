#TODO vectorize constraints and links
mutable struct Combined3{T,N,Nc,Nl,C1,C2,C3,L1,L2,L3} <: Constraint{T,N,Nc,Nl}
    id::Int64
    linkids::SVector{Nl,Int64}

    constr1::C1
    constr2::C2
    constr3::C3
    link1::L1
    link2::L2
    link3::L3

    s0::SVector{N,T}
    s1::SVector{N,T}

    function Combined3(c1, c2, c3)
        constr1,l1,l2 = c1
        constr2,l3,l4 = c2
        constr3,l5,l6 = c3
        constr = [constr1;constr2;constr3]
        links = unique([l1;l2;l3;l4;l5;l6])

        T = constr1.T
        Nl = length(links)
        Nc = constr1.Nc + constr2.Nc + constr3.Nc

        id = getGlobalID()
        linkids = unique([links[i].id for i=1:Nl])

        s0 = @SVector zeros(T,Nc)
        s1 = @SVector zeros(T,Nc)

        new{T,Nc,3,Nl,typeof(constr[1]),typeof(constr[2]),typeof(constr[3]),typeof(links[1]),typeof(links[2]),typeof(links[3])}(id,linkids,constr...,links...,s0,s1)
    end
end

g(C::Combined3) = [g(C.constr1,C.link1,C.link2);g(C.constr2,C.link1,C.link2);g(C.constr3,C.link1,C.link3)]

function ∂g∂pos(C::Combined3,L::Link)
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

function ∂g∂vel(C::Combined3,L::Link)
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
