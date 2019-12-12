abstract type Constraint{T,Nc,Nc²,Nl} <: Node{T,Nc} end

struct Combined2{T,Nc,Nc²,Nl,C1,C2,L1,L2} <: Constraint{T,Nc,Nc²,Nl}
    constr1::C1
    constr2::C2
    data::NodeData{T,Nc,Nc²}
    link1::L1
    link2::L2

    # linkids::SVector{Nl,Int64}



    function Combined2(c1, c2)
        constr1,link1,link2 = c1
        constr2,link3,link4 = c2
        ids = unique([link1.data.id;link2.data.id;link3.data.id;link4.data.id;])
        links = unique([link1;link2;link3;link4])

        T = constr1.T
        Nc = constr1.Nc+constr2.Nc
        Nc² = Nc^2
        Nl = length(ids)
        data = NodeData{T,Nc}()

        # link1 = constr1.link1
        # link2 = constr1.link2


        new{T,Nc,Nc²,Nl,typeof(constr1),typeof(constr2),typeof(link1),typeof(link2)}(constr1,constr2,data,links...)
    end
end



function Combined2(joint::XMLJoint, link1::Link, link2::Link)
    Combined2(Socket(link1,link2,joint.pids...),Axis(link1,link2,joint.axis))
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, C::Constraint)
    summary(io, C); println(io)
    print(io, "\nConnected links: ")
    show(io, mime, C.linkids)
end

Base.show(io::IO, C::Constraint) = summary(io, C)


@inline g(C::Combined2) = [g(C.constr1,C.link1,C.link2);g(C.constr2,C.link1,C.link2)]

@inline function ∂g∂pos(C::Combined2,L::Link)
    if L.data.id == linkids(C)[1]
        return [∂g∂posa(C.constr1,L,C.link2);∂g∂posa(C.constr2,L,C.link2)]
    else
        return [∂g∂posb(C.constr1,C.link1,L);∂g∂posb(C.constr2,C.link1,L)]
    end
end

@inline function ∂g∂vel(C::Combined2,L::Link)
    if L.data.id == linkids(C)[1]
        return [∂g∂vela(C.constr1,L,C.link2);∂g∂vela(C.constr2,L,C.link2)]
    else
        return [∂g∂velb(C.constr1,C.link1,L);∂g∂velb(C.constr2,C.link1,L)]
    end
end

# @inline ∂g∂posa(C::Combined2) = [∂g∂posa(C.constr1);∂g∂posa(C.constr2)]
# @inline ∂g∂posb(C::Combined2) = [∂g∂posb(C.constr1);∂g∂posb(C.constr2)]
# @inline ∂g∂vela(C::Combined2) = [∂g∂vela(C.constr1);∂g∂vela(C.constr2)]
# @inline ∂g∂velb(C::Combined2) = [∂g∂velb(C.constr1);∂g∂velb(C.constr2)]

# @inline Gtλ(L::Link,C::Constraint) = L.data.id==C.linkids[1] ? ∂g∂posa(C)'*C.data.s1 : ∂g∂posb(C)'*C.data.s1
@inline function Gtλ(C::Combined2,L::Link)
    ∂g∂pos(C,L)'*C.data.s1
    # if L.data.id==C.linkids[1]
    #     return ∂g∂pos(C,L,C.link2)'*C.data.s1
    # else
    #     ∂g∂posb(C,C.link2,L)'*C.data.s1
    # end
end

linkids(C::Combined2) = SVector{2,Int64}(C.link1.data.id,C.link2.data.id)
