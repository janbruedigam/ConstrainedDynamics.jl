@inline function setD!(link::Link{T}) where T
    dynT,dynR = ∂dyn∂vel(link)
    Z = @SMatrix zeros(T,3,3)

    link.data.D = [[dynT Z];[Z dynR]]
    return nothing
end

@inline setD!(C::Constraint{T,Nc}) where {T,Nc} = (C.data.D = @SMatrix zeros(T,Nc,Nc); nothing)
@inline function updateD!(node::Node,child::Node,fillin::FillIn)
    node.data.D -= fillin.data.JL*child.data.D*fillin.data.JU
    return nothing
end
@inline invertD!(node) = (d = node.data; d.Dinv = inv(d.D); nothing)

@inline function setJ!(L::Link,C::Constraint,F::FillIn)
    data = F.data

    if L.data.id==C.linkids[1]
        data.JL = ∂g∂vela(C)
        data.JU = -∂g∂posa(C)'
    else
        data.JL = ∂g∂velb(C)
        data.JU = -∂g∂posb(C)'
    end

    return nothing
end

@inline function setJ!(C::Constraint,L::Link,F::FillIn)
    data = F.data

    if L.data.id==C.linkids[1]
        data.JL = -∂g∂posa(C)'
        data.JU = ∂g∂vela(C)
    else
        data.JL = -∂g∂posb(C)'
        data.JU = ∂g∂velb(C)
    end

    return nothing
end

@inline function updateJ!(node::Node,F::FillIn)
    d = F.data
    d.JL = d.JL*node.data.Dinv
    d.JU = node.data.Dinv*d.JU
    return nothing
end

@inline setSol!(link::Link) = (link.data.ŝ = dynamics(link); nothing)
@inline setSol!(C::Constraint) = (C.data.ŝ = g(C); nothing)

# (A) For extended equations
# @inline addGtλ!(L::Link,C::Constraint) = (L.data.ŝ -= Gtλ(L,C); nothing)
@inline addλ0!(C::Constraint) = (C.data.ŝ += C.data.s0; nothing)

@inline function LSol!(node::Node,child::Node,fillin::FillIn)
    node.data.ŝ -= fillin.data.JL*child.data.ŝ
    return nothing
end
@inline DSol!(node) = (d = node.data; d.ŝ = d.Dinv*d.ŝ; nothing)
@inline function USol!(node::Node,parent::Node,fillin::FillIn)
    node.data.ŝ -= fillin.data.JU*parent.data.ŝ
    return nothing
end


@inline function setNormf!(link::Link,robot::Robot)
    data = link.data
    data.f = dynamics(link)
    for (cid,isconnected) in enumconnected(robot.graph,data.id)
        isconnected && GtλTof!(link,getnode(robot,cid))
    end
    data.normf = data.f'*data.f
    return nothing
end

@inline function setNormf!(C::Constraint,robot::Robot)
    data = C.data
    data.f = g(C)
    data.normf = data.f'*data.f
    return nothing
end

@inline GtλTof!(L::Link,C::Constraint) = (L.data.f -= Gtλ(L,C); nothing)
