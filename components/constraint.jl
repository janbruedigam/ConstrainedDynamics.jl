abstract type Constraint{T,Nc,N,Nc²,NcN,Nl} <: Node{T,Nc} end

struct Combined{T,Nc,N,Nc²,NcN,Nl,C1,C2} <: Constraint{T,Nc,N,Nc²,NcN,Nl}
    constr1::C1
    constr2::C2

    linkids::SVector{Nl,Int64}

    data::NodeData{T,Nc,N,Nc²,NcN}

    function Combined(constr1::C1, constr2::C2) where {C1,C2}
        ids = unique([linkids(constr1);linkids(constr2)])

        T = constr1.T
        Nc = constr1.Nc+constr2.Nc
        N = 6
        Nc² = Nc^2
        NcN = Nc*N
        Nl = length(ids)
        data = NodeData{T,Nc,N}()

        new{T,Nc,N,Nc²,NcN,Nl,C1,C2}(constr1,constr2,ids,data)
    end
end



function Combined(joint::XMLJoint, link1::Link, link2::Link)
    Combined(Socket(link1,link2,joint.pids...),Axis(link1,link2,joint.axis))
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, C::Constraint)
    summary(io, C); println(io)
    print(io, "\nConnected links: ")
    show(io, mime, linkids(C))
end

Base.show(io::IO, C::Constraint) = summary(io, C)


@inline g(C::Combined) = [g(C.constr1);g(C.constr2)]

@inline ∂g∂posa(C::Combined) = [∂g∂posa(C.constr1);∂g∂posa(C.constr2)]
@inline ∂g∂posb(C::Combined) = [∂g∂posb(C.constr1);∂g∂posb(C.constr2)]
@inline ∂g∂vela(C::Combined) = [∂g∂vela(C.constr1);∂g∂vela(C.constr2)]
@inline ∂g∂velb(C::Combined) = [∂g∂velb(C.constr1);∂g∂velb(C.constr2)]

#TODO metaprogramming
# @inline linkids(C::Combined{T,Nc,N,Nc²,NcN,1}) where {T,Nc,N,Nc²,NcN} = @SVector [C.constr1.link.data.id]
# @inline linkids(C::Combined{T,Nc,N,Nc²,NcN,2}) where {T,Nc,N,Nc²,NcN} = @SVector [C.constr1.link1.data.id, C.constr1.link2.data.id]

#TODO proper!
function getlinks(C::Combined)
    [C.constr1.link1;C.constr1.link1]
end

@inline Gtλ(L::Link,C::Constraint) = L.data.id==C.linkids[1] ? ∂g∂posa(C)'*C.data.s1 : ∂g∂posb(C)'*C.data.s1
