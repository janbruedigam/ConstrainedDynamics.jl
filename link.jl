using StaticArrays

mutable struct Link{T,No,Np}
    # Properties Assigned when added to a robot
    id::MVector{1,Int64}
    dt::MVector{1,T}
    g::MVector{1,T}

    # Inertia Properties
    m::T
    J::SMatrix{3,3,T}

    # Position, quaternion, force, torque
    x::MVector{No,SVector{3,T}}
    q::MVector{No,Quaternion{T}}
    F::MVector{No,SVector{3,T}}
    τ::MVector{No,SVector{3,T}}

    # New velocity, new angular velocity
    vnew::MVector{3,T}
    ωnew::MVector{3,T}

    # Attachment points for joints/constraints
    p::SVector{Np,SVector{3,T}}

    # Dynamics
    dynT::Function
    dynR::Function

    # Gradients
    ∂dynT∂vnew::Function
    ∂dynR∂ωnew::Function
end

function Base.show(io::IO, link::Link{T,No,Np}) where {T,No,Np}
    heading = string("Link{",T,",",No,",",Np,"}:")
    id = string("\n Unique ID (id): ",link.id[1])
    m = string("\n Mass (m): ",link.m)
    J = string("\n Intertia (J): ",link.J)
    p = string("\n Attachment points (p): ",link.p...)

    print(io,heading,id,m,J,p)
end

function Link{T,No}(m::T, J::Array{T,2}, p::Array{Array{T,1},1}) where {T,No}
    Np = length(p)

    J = convert(SMatrix{3,3,T}, J)
    p = convert(SVector{Np,SVector{3,T}},p)

    x = mrepeat(sbracket(SVector{3,T}(0,0,0)),No)
    q = mrepeat(sbracket(Quaternion{T}()),No)
    F = mrepeat(sbracket(SVector{3,T}(0,0,0)),No)
    τ = mrepeat(sbracket(SVector{3,T}(0,0,0)),No)

    vnew = MVector{3,T}(0,0,0)
    ωnew = MVector{3,T}(0,0,0)

    id = MVector{1,Int64}(0)
    dt = MVector{1,T}(0)
    g = MVector{1,T}(0)

    dynT = function()
        ez = SVector{3,T}(0,0,1)
        v1 = (x[2]-x[1])./dt[1]

        return m.*((vnew - v1)./dt[1] + g[1].*ez) - F[2]
    end
    dynR = function()
        ω1 = 2/dt[1].*Vmat(T)*LTmat(q[1])*q[2]
        sq1 = sqrt(4/dt[1]^2 - ω1'*ω1)*J*ω1
        sq2 = sqrt(4/dt[1]^2 - ωnew'*ωnew)*J*ωnew
        return sq2 + skew(ωnew)*(J*ωnew) - (sq1 - skew(ω1)*(J*ω1)) - 2*τ[2]
    end

    ∂dynT∂vnew = () -> m/dt[1].*SMatrix{3,3,T}(I)
    ∂dynR∂ωnew = () -> sqrt(4/dt[1]^2 - ωnew'*ωnew)*J - J*ωnew*ωnew'/sqrt(4/dt[1]^2 - ωnew'*ωnew) + skew(ωnew)*J - skew(J*ωnew)

    Link{T,No,Np}(id,dt,g,m,J,x,q,F,τ,vnew,ωnew,p,dynT,dynR,∂dynT∂vnew,∂dynR∂ωnew)
end
function Link{T}(m::T, J::Array{T,2}, p::Array{Array{T,1},1}) where T
    No = 2
    Link{T,No}(m, J, p)
end
Link(m::T, J::Array{T,2}, p::Array{Array{T,1},1}) where T = Link{T}(m, J, p)
Link(m::T, J::Array{T,2}) where T = Link(m, J, [zeros(T,3)])
Link(T::Type) = Link(one(T), Matrix{T}(I,3,3), [zeros(T,3)])


function setInit!(link::Link{T,No}; x::AbstractVector{T}=zeros(T,3), q::Quaternion{T}=Quaternion{T}(), F::AbstractVector{T}=zeros(T,3), τ::AbstractVector{T}=zeros(T,3)) where {T,No}
    for i=1:No
        link.x[i] = convert(SVector{3,T},x)
        link.q[i] = q
        link.F[i] = convert(SVector{3,T},F)
        link.τ[i] = convert(SVector{3,T},τ)
    end
end

function initialPosition(link1::Link{T},link2::Link{T},q2::Quaternion{T},pId::AbstractVector{Int64},k::Int) where T
    p1 = Quaternion(link1.p[pId[1]])
    p2 = Quaternion(link2.p[pId[2]])

    x2 = link1.x[k] + rotate(p1,link1.q[k]) - rotate(p2,q2)
end

getv1(link) = (link.x[2]-link.x[1])./link.dt[1]
getω1(link::Link{T}) where T = 2/link.dt[1].*Vmat(T)*LTmat(link.q[1])*link.q[2]
getx3(link) = link.vnew*link.dt[1] + link.x[2]
getq3(link) = Quaternion(link.dt[1]/2 .*Lmat(link.q[2])*ωbar(link))
derivωbar(link::Link{T}) where T = [-(link.ωnew./(ωbar(link)[1]))';SMatrix{3,3,T}(I)]
ωbar(link::Link) = Quaternion(sqrt(4/link.dt[1]^2 - link.ωnew'*link.ωnew),link.ωnew)
ωbar(ω::AbstractVector{T},dt::T) where T = Quaternion(sqrt(4/dt^2 - ω'*ω),ω)
