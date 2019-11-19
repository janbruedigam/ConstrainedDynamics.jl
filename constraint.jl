using StaticArrays

mutable struct Constraint{T,Nc}
    type::String

    λ::SVector{Nc,T}

    linkids::Function

    g::Function

    ∂g∂xa::Function
    ∂g∂xb::Function
    ∂g∂qa::Function
    ∂g∂qb::Function

    ∂g∂va::Function
    ∂g∂vb::Function
    ∂g∂ωa::Function
    ∂g∂ωb::Function
end

function Base.show(io::IO, constraint::Constraint{T,Nc}) where {T,Nc}
    heading = string(constraint.type,"-constraint{",T,",",Nc,"} (",Nc," constraints):")
    links = string("\n Connected links (linkids): ",constraint.linkids())

    print(io,heading,links)
end


function SocketConstraint(link1::Link{T,No},link2::Link{T,No},pId::AbstractVector{Int64}) where {T,No}
    Nc = 3
    type = "socket(2)"
    λ = @SVector zeros(T,Nc)

    linkids = () -> (link1.id[1],link2.id[1])

    p1 = Quaternion(link1.p[pId[1]])
    p2 = Quaternion(link2.p[pId[2]])

    g = () -> getx3(link1) + rotate(p1,getq3(link1)) - (getx3(link2) + rotate(p2,getq3(link2)))

    ∂g∂xa = () -> SMatrix{3,3,T}(I)
    ∂g∂xb = () -> -SMatrix{3,3,T}(I)
    ∂g∂qa = () -> 2*Vmat(T)*RTmat(link1.q[No])*Rmat(p1)*Lmat(link1.q[No])*VTmat(T)
    ∂g∂qb = () -> -2*Vmat(T)*RTmat(link2.q[No])*Rmat(p2)*Lmat(link2.q[No])*VTmat(T)

    ∂g∂va = () -> link1.dt[1]*SMatrix{3,3,T}(I)
    ∂g∂vb = () -> -link2.dt[1]*SMatrix{3,3,T}(I)
    ∂g∂ωa = () -> 2*link1.dt[1]^2/4 .*Vmat(T)*RTmat(link1.q[No])*Lmat(link1.q[No])*RTmat(ωbar(link1))*Rmat(p1)*derivωbar(link1)
    ∂g∂ωb = () -> -2*link2.dt[1]^2/4 .*Vmat(T)*RTmat(link2.q[No])*Lmat(link2.q[No])*RTmat(ωbar(link2))*Rmat(p2)*derivωbar(link2)

    Constraint{T,Nc}(type,λ,linkids,g,∂g∂xa,∂g∂xb,∂g∂qa,∂g∂qb,∂g∂va,∂g∂vb,∂g∂ωa,∂g∂ωb)
end

function SocketConstraint(link::Link{T,No},pId::Int64) where {T,No}
    Nc = 3
    type = "socket(1)"
    λ = @SVector zeros(T,Nc)

    linkids = () -> (link.id[1], link.id[1])

    p = Quaternion(link.p[pId])

    g = () -> getx3(link) + rotate(p,getq3(link))

    ∂g∂xa = () -> SMatrix{3,3,T}(I)
    ∂g∂xb = () -> @SMatrix zeros(T,3,3)
    ∂g∂qa = () -> 2*Vmat(T)*RTmat(link.q[No])*Rmat(p)*Lmat(link.q[No])*VTmat(T)
    ∂g∂qb = () -> @SMatrix zeros(T,3,3)

    ∂g∂va = () -> link.dt[1]*SMatrix{3,3,T}(I)
    ∂g∂vb = () -> @SMatrix zeros(T,3,3)
    ∂g∂ωa = () -> 2*link.dt[1]^2/4 .*Vmat(T)*RTmat(link.q[No])*Lmat(link.q[No])*RTmat(ωbar(link))*Rmat(p)*derivωbar(link)
    ∂g∂ωb = () -> @SMatrix zeros(T,3,3)

    Constraint{T,Nc}(type,λ,linkids,g,∂g∂xa,∂g∂xb,∂g∂qa,∂g∂qb,∂g∂va,∂g∂vb,∂g∂ωa,∂g∂ωb)
end

function RotationConstraint(link1::Link{T,No},link2::Link{T,No},axis::AbstractVector{T}) where {T,No}
    Nc = 2
    type = "axis"
    λ = @SVector zeros(T,Nc)

    linkids = () -> (link1.id[1],link2.id[1])

    axis = convert(SVector{3,T},axis)
    V12 = (@SMatrix [1 0 0; 0 1 0])*svd(skew(axis)).Vt

    g = () -> V12*Vmat(T)*LTmat(getq3(link1))*getq3(link2)
    # g = () -> skew(axis)*Vmat(T)*LTmat(getq3(link1))*getq3(link2)

    ∂g∂xa = () -> @SMatrix zeros(T,2,3)
    ∂g∂xb = () -> @SMatrix zeros(T,2,3)
    ∂g∂qa = () -> -V12*Vmat(T)*Rmat(link2.q[No])*RTmat(link1.q[No])*VTmat(T)
    ∂g∂qb = () -> V12*Vmat(T)*LTmat(link1.q[No])*Lmat(link2.q[No])*VTmat(T)

    ∂g∂va = () -> @SMatrix zeros(T,2,3)
    ∂g∂vb = () -> @SMatrix zeros(T,2,3)
    ∂g∂ωa = () -> link1.dt[1]^2/4*V12*Vmat(T)*Rmat(ωbar(link2))*Rmat(link2.q[No])*RTmat(link1.q[No])*Tmat(T)*derivωbar(link1)
    ∂g∂ωb = () -> link2.dt[1]^2/4*V12*Vmat(T)*LTmat(ωbar(link1))*LTmat(link1.q[No])*Lmat(link2.q[No])*derivωbar(link2)

    Constraint{T,Nc}(type,λ,linkids,g,∂g∂xa,∂g∂xb,∂g∂qa,∂g∂qb,∂g∂va,∂g∂vb,∂g∂ωa,∂g∂ωb)
end

function RotationConstraint(link::Link{T,No}) where {T,No}
    Nc = 3
    k=2
    type = "rotation-lock"
    λ = @SVector zeros(T,Nc)

    linkids = () ->  (link.id[1], link.id[1])

    g = () -> Vmat(T)*getq3(link)

    ∂g∂xa = () -> @SMatrix zeros(T,3,3)
    ∂g∂xb = () -> @SMatrix zeros(T,3,3)
    ∂g∂qa = () -> Vmat(T)*Lmat(link.q[No])*VTmat(T)
    ∂g∂qb = () -> @SMatrix zeros(T,3,3)

    ∂g∂va = () -> @SMatrix zeros(T,3,3)
    ∂g∂vb = () -> @SMatrix zeros(T,3,3)
    ∂g∂ωa = () -> link.dt[1]/2 .*Vmat(T)*Lmat(link.q[No])*derivωbar(link)
    ∂g∂ωb = () -> @SMatrix zeros(T,3,3)

    Constraint{T,Nc}(type,λ,linkids,g,∂g∂xa,∂g∂xb,∂g∂qa,∂g∂qb,∂g∂va,∂g∂vb,∂g∂ωa,∂g∂ωb)
end
