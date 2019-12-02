include("../node.jl")
include("../link.jl")

struct Axis{T,Nc,N,Nc²,NcN,Nl,L1,L2} <: Constraint{T,Nc,N,Nc²,NcN,Nl}
    link1::L1
    link2::L2

    V12::SMatrix{2,3,T,6}

    data::NodeData{T,Nc,N,Nc²,NcN}
end

function Axis(link1::Link{T,N},link2::Link{T,N},axis::AbstractVector{T}) where {T,N}
    Nc = 2
    Nc² = Nc^2
    NcN = Nc*N
    Nl = 2
    V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt
    data = NodeData{T,Nc,N}()

    Axis{T,Nc,N,Nc²,NcN,Nl,typeof(link1),typeof(link2)}(link1,link2,V12,data)
end


@inline g(C::Axis) = C.V12*Vmat(LTmat(getq3(C.link1))*getq3(C.link2))

@inline function ∂g∂posa(C::Axis{T}) where T
    X = @SMatrix zeros(T,2,3)

    link1 = C.link1
    link2 = C.link2
    No = link1.No
    R = -C.V12*Vmat(VTmat(Rmat(link2.q[No])*RTmat(link1.q[No])))

    return [X R]
end

@inline function ∂g∂posb(C::Axis{T}) where T
    X = @SMatrix zeros(T,2,3)

    link1 = C.link1
    link2 = C.link2
    No = link1.No
    R = C.V12*Vmat(VTmat(LTmat(link1.q[No])*Lmat(link2.q[No])))

    return [X R]
end

@inline function ∂g∂vela(C::Axis{T}) where T
    V = @SMatrix zeros(T,2,3)

    link1 = C.link1
    link2 = C.link2
    No = link1.No
    Ω = link1.dt^2/4*C.V12*Vmat(Rmat(ωbar(link2))*Rmat(link2.q[No])*RTmat(link1.q[No])*Tmat(T))*derivωbar(link1)

    return [V Ω]
end

@inline function ∂g∂velb(C::Axis{T}) where T
    V = @SMatrix zeros(T,2,3)

    link1 = C.link1
    link2 = C.link2
    No = link1.No
    Ω = link2.dt^2/4*C.V12*Vmat(LTmat(ωbar(link1))*LTmat(link1.q[No])*Lmat(link2.q[No]))*derivωbar(link2)

    return [V Ω]
end
