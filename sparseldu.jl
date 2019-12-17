@inline function setD!(diagonal,link::Link{T}) where T
    # μ = 1e-4
    diagonal.D = ∂dyn∂vel(link) #+ SMatrix{6,6,Float64,36}(μ*I)
    return nothing
end

@inline function setD!(diagonal,C::Constraint{T,Nc}) where {T,Nc}
    diagonal.D = @SMatrix zeros(T,Nc,Nc)
    return nothing
end


# TODO pass in the two connected links
@inline function setJ!(F::OffDiagonalEntry,C::Constraint,L::Link)
    # data = F.data

    F.JL = ∂g∂vel(C,L)
    F.JU = -∂g∂pos(C,L)'

    return nothing
end

@inline function setJ!(F::OffDiagonalEntry,L::Link,C::Constraint)
    F.JL = -∂g∂pos(C,L)'
    F.JU = ∂g∂vel(C,L)

    return nothing
end

@inline function setJ!(F::OffDiagonalEntry{T,N1,N2},C1::Constraint,C2::Constraint) where {T,N1,N2}
    F.JL = @SMatrix zeros(T,N2,N1)
    F.JU = @SMatrix zeros(T,N1,N2)

    return nothing
end

@inline function setJ!(entry::OffDiagonalEntry{T,N1,N2}) where {T,N1,N2}
    entry.JL = @SMatrix zeros(T,N2,N1)
    entry.JU = @SMatrix zeros(T,N1,N2)
    return nothing
end

@inline setSol!(diagonal,link::Link) = (diagonal.ŝ = +(dynamics(link)); nothing)
@inline setSol!(diagonal,C::Constraint) = (diagonal.ŝ = g(C); nothing)

# (A) For extended equations
# @inline addGtλ!(L::Link,C::Constraint) = (L.data.ŝ -= Gtλ(L,C); nothing)
@inline addλ0!(diagonal,C::Constraint) = (diagonal.ŝ += C.data.s0; nothing)


@inline function setNormf!(link::Link,robot::Robot)
    data = link.data
    graph = robot.graph
    data.f = dynamics(link)
    for (cind,isconnected) in enumerate(graph.adjacency[graph.dict[link.id]])
        isconnected && GtλTof!(robot.nodes[robot.dict[graph.rdict[cind]]],link)
    end
    data.normf = data.f'*data.f
    return nothing
end

@inline function setNormf!(C::Constraint,robot::Robot)
    data = C.data
    data.f = g(C)
    # for i=1:Nl
    #     gf(C.constr[i],C.link1,C.link2,C.datavec[i])
    # end
    data.normf = data.f'*data.f
    return nothing
end

@inline GtλTof!(C::Constraint,L::Link) = (L.data.f -= ∂g∂pos(C,L)'*C.data.s1; nothing)
