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
function setJ!(robot,F::OffDiagonalEntry,C::Constraint,L::Link)
    F.JL = ∂g∂vel(robot,C,L.id)
    F.JU = -∂g∂pos(robot,C,L.id)'
    return
end

function setJ!(robot,F::OffDiagonalEntry,L::Link,C::Constraint)
    F.JL = -∂g∂pos(robot,C,L.id)'
    F.JU = ∂g∂vel(robot,C,L.id)
    return
end

function setJ!(robot,F::OffDiagonalEntry{T,N1,N2},C1::Constraint,C2::Constraint) where {T,N1,N2}
    F.JL = @SMatrix zeros(T,N2,N1)
    F.JU = @SMatrix zeros(T,N1,N2)
    return
end

# function setJ!(robot,entry::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
#     entry.JL = @SMatrix zeros(T,N2,N1)
#     entry.JU = @SMatrix zeros(T,N1,N2)
#     return
# end

setSol!(diagonal,link::Link,robot) = (diagonal.ŝ = dynamics(link); nothing)
setSol!(diagonal,C::Constraint,robot) = (diagonal.ŝ = g(robot,C); nothing)


# (A) For extended equations
# addGtλ!(L::Link,C::Constraint) = (L.data.ŝ -= Gtλ(L,C); nothing)
addλ0!(diagonal,C::Constraint) = (diagonal.ŝ += C.s0; nothing)


function normf(link::Link,robot::Robot)
    graph = robot.graph
    id = link.id
    link.f = dynamics(link)

    for cid in connections(graph,id)
        cid == -1 && break
        GtλTof!(robot,getnode(robot,cid),link)
    end
    f = link.f
    return dot(f,f)
end

function normf(C::Constraint,robot::Robot)
    f = g(robot,C)
    return dot(f,f)
end

GtλTof!(robot,C::Constraint,L::Link) = (L.f -= ∂g∂pos(robot,C,L.id)'*C.s1; nothing)
