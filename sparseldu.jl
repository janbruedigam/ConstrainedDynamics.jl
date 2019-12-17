function setD!(diagonal,link::Link{T}) where T
    # μ = 1e-4
    diagonal.D = ∂dyn∂vel(link) #+ SMatrix{6,6,Float64,36}(μ*I)
    return
end

function setD!(diagonal::DiagonalEntry{T,N},::Constraint) where {T,N}
    diagonal.D = @SMatrix zeros(T,N,N)
    return
end


# TODO pass in the two connected links
function setJ!(F::OffDiagonalEntry,C::Constraint,L::Link)
    F.JL = ∂g∂vel(C,L)
    F.JU = -∂g∂pos(C,L)'
    return
end

function setJ!(F::OffDiagonalEntry,L::Link,C::Constraint)
    F.JL = -∂g∂pos(C,L)'
    F.JU = ∂g∂vel(C,L)
    return
end

function setJ!(F::OffDiagonalEntry{T,N1,N2},C1::Constraint,C2::Constraint) where {T,N1,N2}
    F.JL = @SMatrix zeros(T,N2,N1)
    F.JU = @SMatrix zeros(T,N1,N2)
    return
end

function setJ!(entry::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
    entry.JL = @SMatrix zeros(T,N2,N1)
    entry.JU = @SMatrix zeros(T,N1,N2)
    return
end

setSol!(diagonal,link::Link) = (diagonal.ŝ = dynamics(link); nothing)
setSol!(diagonal,C::Constraint) = (diagonal.ŝ = g(C); nothing)


# (A) For extended equations
# addGtλ!(L::Link,C::Constraint) = (L.data.ŝ -= Gtλ(L,C); nothing)
addλ0!(diagonal,C::Constraint) = (diagonal.ŝ += C.s0; nothing)


function normf(link::Link,robot::Robot)
    graph = robot.graph
    link.f = dynamics(link)
    for cid in connections(graph,link.id)
        cid == -1 && break
        GtλTof!(getnode(robot,cid),link)
    end
    f = link.f
    return dot(f,f)
end

function normf(C::Constraint,robot::Robot)
    f = g(C)
    return dot(f,f)
end

GtλTof!(C::Constraint,L::Link) = (L.f -= ∂g∂pos(C,L)'*C.s1; nothing)
