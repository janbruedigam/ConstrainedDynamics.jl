struct Axis{T,Nc,Nl,L1,L2} <: Joint{T,Nc,Nl}
    link1::L1
    link2::L2

    V12::SMatrix{2,3,T,6}

    function Axis(link1::Link{T},link2::Link{T},axis::AbstractVector{T}) where T
        Nc = 2
        Nl = 2
        V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt

        new{T,Nc,Nl,typeof(link1),typeof(link2)}(link1,link2,V12)
    end
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
