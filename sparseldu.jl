# @inline function setD!(link::Link{T}) where T
#     μ = 1e-4
#     dynT,dynR = ∂dyn∂vel(link)
#     Z = @SMatrix zeros(T,3,3)
#
#     link.data.D = [[dynT Z];[Z dynR]] #+ SMatrix{6,6,Float64,36}(μ*I)
#     return nothing
# end
#
# @inline setD!(C::Constraint{T,Nc}) where {T,Nc} = (C.data.D = @SMatrix zeros(T,Nc,Nc); nothing)

@inline function setD!(diagonal,link::Link{T}) where T
    μ = 1e-4
    dynT,dynR = ∂dyn∂vel(link)
    Z = @SMatrix zeros(T,3,3)

    diagonal.D = [[dynT Z];[Z dynR]] #+ SMatrix{6,6,Float64,36}(μ*I)
    return nothing
end

@inline function setD!(diagonal,C::Constraint{T,Nc}) where {T,Nc}
    diagonal.D = @SMatrix zeros(T,Nc,Nc)
    return nothing
end


@inline function updateD!(diagonal::DiagonalEntry,child::DiagonalEntry,fillin::OffDiagonalEntry)
    diagonal.D -= fillin.JL*child.D*fillin.JU
    return nothing
end
@inline invertD!(diagonal) = (diagonal.Dinv = +inv(diagonal.D); nothing)

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

@inline function updateJ1!(F::OffDiagonalEntry,diagonal::DiagonalEntry,gcfillin::OffDiagonalEntry,cgcfillin::OffDiagonalEntry)
    F.JL -= gcfillin.JL*diagonal.D*cgcfillin.JU
    F.JU -= cgcfillin.JL*diagonal.D*gcfillin.JU
    return nothing
end

@inline function updateJ2!(F::OffDiagonalEntry,diagonal::DiagonalEntry)
    F.JL = F.JL*diagonal.Dinv
    F.JU = diagonal.Dinv*F.JU
    return nothing
end

@inline setSol!(diagonal,link::Link) = (diagonal.ŝ = diagonal.ŝ -diagonal.ŝ +dynamics(link); nothing)
@inline setSol!(diagonal,C::Constraint) = (diagonal.ŝ = g(C); nothing)

# (A) For extended equations
# @inline addGtλ!(L::Link,C::Constraint) = (L.data.ŝ -= Gtλ(L,C); nothing)
@inline addλ0!(diagonal,C::Constraint) = (diagonal.ŝ += C.data.s0; nothing)

@inline function LSol!(diagonal::DiagonalEntry,child::DiagonalEntry,fillin::OffDiagonalEntry)
    diagonal.ŝ -= fillin.JL*child.ŝ
    return nothing
end

#TODO whats going on here ???
@inline DSol!(diagonal) = (diagonal.ŝ = diagonal.ŝ - diagonal.ŝ + diagonal.Dinv*diagonal.ŝ; nothing)
@inline function USol!(diagonal::DiagonalEntry,parent::DiagonalEntry,fillin::OffDiagonalEntry)
    diagonal.ŝ -= fillin.JU*parent.ŝ
    return nothing
end


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
