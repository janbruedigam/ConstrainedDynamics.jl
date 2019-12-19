@inline function setDandŝ!(diagonal,link::Link,robot)
    diagonal.D = ∂dyn∂vel(robot, link)
    diagonal.ŝ = dynamics(robot, link)
    return
end

@inline function setDandŝ!(diagonal::DiagonalEntry{T,N},C::Constraint,robot) where {T,N}
    diagonal.D = @SMatrix zeros(T,N,N)
    diagonal.ŝ = g(robot,C)
    return
end

# TODO pass in the two connected links
@inline function setJ!(robot,F::OffDiagonalEntry,linkid::Int64,C::Constraint)
    F.JL = -∂g∂pos(robot,C,linkid)'
    F.JU = ∂g∂vel(robot,C,linkid)
    return
end

@inline function setJ!(robot,F::OffDiagonalEntry,C::Constraint,linkid::Int64)
    F.JL = ∂g∂vel(robot,C,linkid)
    F.JU = -∂g∂pos(robot,C,linkid)'
    return
end

@inline function setJ!(F::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
    F.JL = @SMatrix zeros(T,N2,N1)
    F.JU = F.JL'
    return
end

addλ0!(diagonal,C::Constraint) = (diagonal.ŝ += C.s0; return)

# function normf(link::Link{T},robot::Robot) where T
#     f = dynamics(robot,link)
#     return dot(f,f)
# end

@inline function normf(link::Link{T},robot::Robot) where T
    graph = robot.graph
    id = link.id
    link.f = dynamics(robot,link)

    for cid in connections(graph,id)
        cid == -1 && break
        GtλTof!(robot,getconstraint(robot,cid),link)
    end

    f = link.f
    return dot(f,f)
end

@inline function normf(C::Constraint,robot::Robot)
    f = g(robot,C)
    return dot(f,f)
end

GtλTof!(robot,C::Constraint,L::Link) = (L.f -= ∂g∂pos(robot,C,L.id)'*C.s1; return)
