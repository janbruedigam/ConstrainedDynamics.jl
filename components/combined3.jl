#TODO vectorize constraints and links
struct Combined3{T,Nc,Nc²,Nl,C1,C2,C3,L1,L2,L3} <: Constraint{T,Nc,Nc²,Nl}
    constr1::C1
    constr2::C2
    constr3::C3
    link1::L1
    link2::L2
    link3::L3
    data::NodeData{T,Nc,Nc²}

    function Combined3(c1, c2, c3)
        constr1,l1,l2 = c1
        constr2,l3,l4 = c2
        constr3,l5,l6 = c3
        constr = [constr1;constr2;constr3]
        links = unique([l1;l2;l3;l4;l5;l6])

        T = constr1.T
        Nc = constr1.Nc+constr2.Nc+constr3.Nc
        Nc² = Nc^2
        Nl = length(links)

        ids = unique([links[i].data.id for i=1:Nl])

        data = NodeData{T,Nc}()

        new{T,Nc,Nc²,Nl,typeof(constr[1]),typeof(constr[2]),typeof(constr[3]),typeof(links[1]),typeof(links[2]),typeof(links[3])}(constr...,links...,data)
    end
end


function Combined3(joint::XMLJoint, link1::Link, link2::Link)
    Combined3(Socket(link1,link2,joint.pids...),Axis(link1,link2,joint.axis))
end

@inline g(C::Combined3) = [g(C.constr1,C.link1,C.link2);g(C.constr2,C.link1,C.link2)]

@inline function ∂g∂pos(C::Combined3,L::Link)
    if L.data.id == linkids(C)[1]
        return [∂g∂posa(C.constr1,L,C.link2);∂g∂posa(C.constr2,L,C.link2)]
    else
        return [∂g∂posb(C.constr1,C.link1,L);∂g∂posb(C.constr2,C.link1,L)]
    end
end

@inline function ∂g∂vel(C::Combined3,L::Link)
    if L.data.id == linkids(C)[1]
        return [∂g∂vela(C.constr1,L,C.link2);∂g∂vela(C.constr2,L,C.link2)]
    else
        return [∂g∂velb(C.constr1,C.link1,L);∂g∂velb(C.constr2,C.link1,L)]
    end
end

linkids(C::Combined3) = SVector{2,Int64}(C.link1.data.id,C.link2.data.id)
# @generated function linkids(C::Combined3{T,Nc,Nc²,Nl}) where {T,Nc,Nc²,Nl}
#     ids = [:(C.links[$i].data.id) for i=1:Nl]
#     :(SVector{Nl,Int64}($(ids...)))
# end
