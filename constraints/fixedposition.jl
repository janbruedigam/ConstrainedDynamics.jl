#TODO update

struct FixedPosition{T,Nc,N,Nc²,NcN,Nl,L} <: Constraint{T,Nc,N,Nc²,NcN,Nl}
    link::L

    pid::Int64

    data::NodeData{T,Nc,N,Nc²,NcN}
end

function FixedPosition(link::Link{T,N},pid) where {T,N}
    Nc = 3
    Nc² = Nc^2
    NcN = Nc*N
    Nl = 1
    data = NodeData{T,Nc,N}()

    FixedPosition{T,Nc,N,Nc²,NcN,Nl,typeof(link)}(link,pid,data)
end


@inline function g(C::FixedPosition)
    link = C.link
    getx3(link) + rotate(link.p[C.pid],getq3(link))
end

@inline function ∂g∂posa(C::FixedPosition{T}) where T
    link = C.link
    X = SMatrix{3,3,T,9}(I)

    q = link.q[link.No]
    R = 2*Vmat(VTmat(RTmat(q)*Rmat(Quaternion(link.p[C.pid]))*Lmat(q)))

    return [X R]
end

@inline function ∂g∂vela(C::FixedPosition{T}) where T
    link = C.link
    V = SMatrix{3,3,T,9}(link.dt*I)

    q = link.q[link.No]
    Ω = 2*link.dt^2/4*Vmat(RTmat(q)*Lmat(q)*RTmat(ωbar(link))*Rmat(Quaternion(link.p[C.pid])))*derivωbar(link)

    return [V Ω]
end
