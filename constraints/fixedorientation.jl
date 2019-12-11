#TODO update

struct FixedOrientation{T,Nc,N,Nc²,NcN,Nl,L} <: Constraint{T,Nc,N,Nc²,NcN,Nl}
    link::L

    data::NodeData{T,Nc,N,Nc²,NcN}
end

function FixedOrientation(link::Link{T,N}) where {T,N}
    Nc = 3
    Nc² = Nc^2
    NcN = Nc*N
    Nl = 1
    data = NodeData{T,Nc,N}()

    FixedOrientation{T,Nc,N,Nc²,NcN,Nl,typeof(link)}(link,data)
end


@inline g(C::FixedOrientation{T}) where T = Vmat(getq3(C.link))

@inline function ∂g∂posa(C::FixedOrientation{T}) where T
    link = C.link
    X = @SMatrix zeros(T,3,3)

    R = Vmat(VTmat(Lmat(link.q[link.No])))

    return [X R]
end

@inline function ∂g∂vela(C::FixedOrientation{T}) where T
    link = C.link
    V = @SMatrix zeros(T,3,3)

    Ω = link.dt/2*Vmat(Lmat(link.q[link.No])*derivωbar(link))

    return [V Ω]
end
