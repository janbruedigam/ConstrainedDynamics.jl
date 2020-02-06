abstract type AbstractLink{T} <: Component{T} end

mutable struct Link{T} <: AbstractLink{T}
    id::Int64

    m::T
    J::SMatrix{3,3,T,9}

    x::Vector{SVector{3,T}}
    q::Vector{Quaternion{T}}

    F::Vector{SVector{3,T}}
    τ::Vector{SVector{3,T}}

    s0::SVector{6,T}
    s1::SVector{6,T}
    f::SVector{6,T}

    function Link(m::T,J::AbstractArray{T,2}) where T
        x = [@SVector zeros(T,3)]
        q = [Quaternion{T}()]

        F = [@SVector zeros(T,3)]
        τ = [@SVector zeros(T,3)]

        s0 = @SVector zeros(T,6)
        s1 = @SVector zeros(T,6)
        f = @SVector zeros(T,6)

        new{T}(getGlobalID(),m,J,x,q,F,τ,s0,s1,f)
    end

    function Link(shape::Shape)
        L = Link(shape.m,shape.J)
        push!(shape.links,L)
        return L
    end
end

mutable struct Origin{T} <: AbstractLink{T}
    id::Int64

    Origin{T}() where T = new{T}(getGlobalID())
end


Base.length(C::Link) = 6

function setInit!(link::Link{T};
    x::AbstractVector=zeros(T,3),
    q::Quaternion=Quaternion{T}(),
    F::AbstractVector=zeros(T,3),
    τ::AbstractVector=zeros(T,3)) where T

    link.x[1] = convert(SVector{3,T},x)
    link.q[1] = q
    link.F[1] = convert(SVector{3,T},F)
    link.τ[1] = convert(SVector{3,T},τ)
    return
end

function setInit!(link1::Link{T}, link2::Link{T}, p1::AbstractVector, p2::AbstractVector;
    q::Quaternion=Quaternion{T}(),
    F::AbstractVector=zeros(T,3),
    τ::AbstractVector=zeros(T,3)) where T

    x2 = link1.x[1] + vrotate(p1,link1.q[1]) - vrotate(p2,q)
    setInit!(link2; x=x2, q=q, F=F, τ=τ)
end

function setInit!(link1::Origin{T}, link2::Link{T}, p1::AbstractVector, p2::AbstractVector;
    q::Quaternion=Quaternion{T}(),
    F::AbstractVector=zeros(T,3),
    τ::AbstractVector=zeros(T,3)) where T

    x2 = p1 - vrotate(p2,q)
    setInit!(link2; x=x2, q=q, F=F, τ=τ)
end

@inline getx3(link::Link, dt) = getvnew(link)*dt + link.x[2]
@inline getq3(link::Link, dt) = Quaternion(dt/2*(Lmat(link.q[2])*ωbar(link,dt)))
@inline getv1(link::Link, dt) = (link.x[2]-link.x[1])/dt
@inline function getω1(link::Link, dt)
    q1 = link.q[1]
    q2 = link.q[2]
    2/dt*(q1.s*imag(q2)-q2.s*imag(q1)-cross(imag(q1),imag(q2)))
    # 2/dt*(VLᵀmat(link.q[1])*link.q[2])
end

#TODO use SOneTo and SUnitRange once they are faster
# @inline getvnew(link::Link) = link.s1[SOneTo(3)]
# @inline getωnew(link::Link) = link.s1[SUnitRange(4,6)]
@inline getvnew(link::Link) = link.s1[SVector(1,2,3)]
@inline getωnew(link::Link) = link.s1[SVector(4,5,6)]

@inline function derivωbar(link::Link{T}, dt) where T
    ωnew = getωnew(link)
    msq = -sqrt(4/dt^2 - dot(ωnew,ωnew))
    [ωnew'/msq; SMatrix{3,3,T,9}(I)]
end
@inline function ωbar(link::Link, dt)
    ωnew = getωnew(link)
    Quaternion(sqrt(4/dt^2 - dot(ωnew,ωnew)),ωnew)
end

@inline function dynamics(link::Link{T}, robot) where T
    No = robot.No
    dt = robot.dt

    ezg = SVector{3,T}(0,0,-robot.g)
    dynT = link.m*((getvnew(link) - getv1(link,dt))/dt + ezg) - link.F[No]

    J = link.J
    ω1 = getω1(link,dt)
    ωnew = getωnew(link)
    sq1 = sqrt(4/dt^2 - ω1'*ω1)
    sq2 = sqrt(4/dt^2 - ωnew'*ωnew)
    dynR = skewplusdiag(ωnew,sq2)*(J*ωnew) - skewplusdiag(ω1,sq1)*(J*ω1) - 2*link.τ[No]

    link.f = [dynT;dynR]

    for cid in connections(robot.graph,link.id)
        GtλTof!(link,getconstraint(robot,cid),robot)
    end

    return link.f
end

@inline function ∂dyn∂vel(link::Link{T}, dt) where T
    J = link.J
    ωnew = getωnew(link)
    sq = sqrt(4/dt^2 - ωnew'*ωnew)

    dynT = SMatrix{3,3,T,9}(link.m/dt*I)
    dynR = skewplusdiag(ωnew,sq)*J - J*ωnew*(ωnew'/sq) - skew(J*ωnew)

    Z = @SMatrix zeros(T,3,3)

    return [[dynT; Z] [Z; dynR]]
end
