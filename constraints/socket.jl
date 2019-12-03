struct Socket{T,Nc,N,Nc²,NcN,Nl,L1,L2} <: Constraint{T,Nc,N,Nc²,NcN,Nl}
    link1::L1
    link2::L2

    pids::SVector{2,Int64}

    data::NodeData{T,Nc,N,Nc²,NcN}
end

function Socket(link1::Link{T,N},link2::Link{T,N},pid1::Int64,pid2::Int64) where {T,N}
    Nc = 3
    Nc² = Nc^2
    NcN = Nc*N
    Nl = 2
    pids = SVector(pid1,pid2)
    data = NodeData{T,Nc,N}()

    Socket{T,Nc,N,Nc²,NcN,Nl,typeof(link1),typeof(link2)}(link1,link2,pids,data)
end


@inline function g(C::Socket)
    link1 = C.link1
    link2 = C.link2
    pids = C.pids
    getx3(link1) + rotate(link1.p[pids[1]],getq3(link1)) - (getx3(link2) + rotate(link2.p[pids[2]],getq3(link2)))
end

@inline function ∂g∂posa(C::Socket{T}) where T
    link = C.link1
    X = SMatrix{3,3,T,9}(I)

    q = link.q[link.No]
    R = 2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link.p[C.pids[1]]))*Lmat(q)))

    return [X R]
end

@inline function ∂g∂posb(C::Socket{T}) where T
    link = C.link2
    X = SMatrix{3,3,T,9}(-I)

    q = link.q[link.No]
    R = -2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link.p[C.pids[2]]))*Lmat(q)))

    return [X R]
end

@inline function ∂g∂vela(C::Socket{T}) where T
    link = C.link1
    V = link.dt*SMatrix{3,3,T,9}(I)

    q = link.q[link.No]
    Ω = 2*link.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link))*Rmat(Quaternion(link.p[C.pids[1]])))*derivωbar(link)

    return [V Ω]
end

@inline function ∂g∂velb(C::Socket{T}) where T
    link = C.link2
    V = link.dt*SMatrix{3,3,T,9}(-I)

    q = link.q[link.No]
    Ω = -2*link.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link))*Rmat(Quaternion(link.p[C.pids[2]])))*derivωbar(link)

    return [V Ω]
end
