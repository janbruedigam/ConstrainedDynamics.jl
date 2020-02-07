abstract type AbstractBody{T} <: Component{T} end

mutable struct Body{T} <: AbstractBody{T}
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

    function Body(m::T,J::AbstractArray{T,2}) where T
        x = [@SVector zeros(T,3)]
        q = [Quaternion{T}()]

        F = [@SVector zeros(T,3)]
        τ = [@SVector zeros(T,3)]

        s0 = @SVector zeros(T,6)
        s1 = @SVector zeros(T,6)
        f = @SVector zeros(T,6)

        new{T}(getGlobalID(),m,J,x,q,F,τ,s0,s1,f)
    end

    function Body(shape::Shape)
        L = Body(shape.m,shape.J)
        push!(shape.bodies,L)
        return L
    end
end

mutable struct Origin{T} <: AbstractBody{T}
    id::Int64

    Origin{T}() where T = new{T}(getGlobalID())
end


Base.length(C::Body) = 6

function setInit!(body::Body{T};
    x::AbstractVector=zeros(T,3),
    q::Quaternion=Quaternion{T}(),
    F::AbstractVector=zeros(T,3),
    τ::AbstractVector=zeros(T,3)) where T

    body.x[1] = convert(SVector{3,T},x)
    body.q[1] = q
    body.F[1] = convert(SVector{3,T},F)
    body.τ[1] = convert(SVector{3,T},τ)
    return
end

function setInit!(body1::Body{T}, body2::Body{T}, p1::AbstractVector, p2::AbstractVector;
    q::Quaternion=Quaternion{T}(),
    F::AbstractVector=zeros(T,3),
    τ::AbstractVector=zeros(T,3)) where T

    x2 = body1.x[1] + vrotate(p1,body1.q[1]) - vrotate(p2,q)
    setInit!(body2; x=x2, q=q, F=F, τ=τ)
end

function setInit!(body1::Origin{T}, body2::Body{T}, p1::AbstractVector, p2::AbstractVector;
    q::Quaternion=Quaternion{T}(),
    F::AbstractVector=zeros(T,3),
    τ::AbstractVector=zeros(T,3)) where T

    x2 = p1 - vrotate(p2,q)
    setInit!(body2; x=x2, q=q, F=F, τ=τ)
end

@inline getx3(body::Body, dt) = getvnew(body)*dt + body.x[2]
@inline getq3(body::Body, dt) = Quaternion(dt/2*(Lmat(body.q[2])*ωbar(body,dt)))
@inline getv1(body::Body, dt) = (body.x[2]-body.x[1])/dt
@inline function getω1(body::Body, dt)
    q1 = body.q[1]
    q2 = body.q[2]
    2/dt*(q1.s*imag(q2)-q2.s*imag(q1)-cross(imag(q1),imag(q2)))
    # 2/dt*(VLᵀmat(body.q[1])*body.q[2])
end

#TODO use SOneTo and SUnitRange once they are faster
# @inline getvnew(body::Body) = body.s1[SOneTo(3)]
# @inline getωnew(body::Body) = body.s1[SUnitRange(4,6)]
@inline getvnew(body::Body) = body.s1[SVector(1,2,3)]
@inline getωnew(body::Body) = body.s1[SVector(4,5,6)]

@inline function derivωbar(body::Body{T}, dt) where T
    ωnew = getωnew(body)
    msq = -sqrt(4/dt^2 - dot(ωnew,ωnew))
    [ωnew'/msq; SMatrix{3,3,T,9}(I)]
end
@inline function ωbar(body::Body, dt)
    ωnew = getωnew(body)
    Quaternion(sqrt(4/dt^2 - dot(ωnew,ωnew)),ωnew)
end

@inline function dynamics(body::Body{T}, robot) where T
    No = robot.No
    dt = robot.dt

    ezg = SVector{3,T}(0,0,-robot.g)
    dynT = body.m*((getvnew(body) - getv1(body,dt))/dt + ezg) - body.F[No]

    J = body.J
    ω1 = getω1(body,dt)
    ωnew = getωnew(body)
    sq1 = sqrt(4/dt^2 - ω1'*ω1)
    sq2 = sqrt(4/dt^2 - ωnew'*ωnew)
    dynR = skewplusdiag(ωnew,sq2)*(J*ωnew) - skewplusdiag(ω1,sq1)*(J*ω1) - 2*body.τ[No]

    body.f = [dynT;dynR]

    for cid in connections(robot.graph,body.id)
        GtλTof!(body,getconstraint(robot,cid),robot)
    end

    return body.f
end

@inline function ∂dyn∂vel(body::Body{T}, dt) where T
    J = body.J
    ωnew = getωnew(body)
    sq = sqrt(4/dt^2 - ωnew'*ωnew)

    dynT = SMatrix{3,3,T,9}(body.m/dt*I)
    dynR = skewplusdiag(ωnew,sq)*J - J*ωnew*(ωnew'/sq) - skew(J*ωnew)

    Z = @SMatrix zeros(T,3,3)

    return [[dynT; Z] [Z; dynR]]
end
