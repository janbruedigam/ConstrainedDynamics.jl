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


    function Body(m::T, J::AbstractArray{T,2}) where T
        x = [zeros(T, 3)]
        q = [Quaternion{T}()]

        F = [zeros(T, 3)]
        τ = [zeros(T, 3)]

        s0 = zeros(T, 6)
        s1 = zeros(T, 6)
        f = zeros(T, 6)

        new{T}(getGlobalID(), m, J, x, q, F, τ, s0, s1, f)
    end

    function Body(shape::Shape)
        body = Body(shape.m, shape.J)
        push!(shape.bodyids, body.id)
        return body
    end

    Body(m::T, J::SMatrix{3,3,T,9}, x::Vector{SVector{3,T}}, q::Vector{Quaternion{T}}, F::Vector{SVector{3,T}}, τ::Vector{SVector{3,T}}, 
        s0::SVector{6,T}, s1::SVector{6,T}, f::SVector{6,T}) where T = new{T}(getGlobalID(),m,J,x,q,F,τ,s0,s1,f)
end

mutable struct Origin{T} <: AbstractBody{T}
    id::Int64

    Origin{T}() where T = new{T}(getGlobalID())
end

function Base.deepcopy(b::Body)
    contents = []
    for i = 2:10
        push!(contents, deepcopy(getfield(b, i)))
    end

    Body(contents...)
end

function Base.deepcopy(b::Body, shape::Shape)
    contents = []
    for i = 2:10
        push!(contents, deepcopy(getfield(b, i)))
    end

    body = Body(contents...)
    @assert (body.m == shape.m && body.J == shape.J)
    push!(shape.bodyids, body.id)
    return body
end

Base.length(::Body) = 6

function setInit!(body::Body{T};
    x::AbstractVector = zeros(T, 3),
    q::Quaternion = Quaternion{T}(),
    F::AbstractVector = zeros(T, 3),
    τ::AbstractVector = zeros(T, 3)) where T

    body.x[1] = convert(SVector{3,T}, x)
    body.q[1] = q
    body.F[1] = convert(SVector{3,T}, F)
    body.τ[1] = convert(SVector{3,T}, τ)
    return
end

function setInit!(body1::Body{T}, body2::Body{T}, p1::AbstractVector, p2::AbstractVector;
    q::Quaternion = Quaternion{T}(),
    F::AbstractVector = zeros(T, 3),
    τ::AbstractVector = zeros(T, 3)) where T

    p1 = convert(SVector{3,T}, p1)
    p2 = convert(SVector{3,T}, p2)

    x2 = body1.x[1] + vrotate(p1, body1.q[1]) - vrotate(p2, q)
    setInit!(body2; x = x2, q = q, F = F, τ = τ)
end

function setInit!(body1::Origin{T}, body2::Body{T}, p1::AbstractVector, p2::AbstractVector;
    q::Quaternion = Quaternion{T}(),
    F::AbstractVector = zeros(T, 3),
    τ::AbstractVector = zeros(T, 3)) where T

    p2 = convert(SVector{3,T}, p2)

    x2 = p1 - vrotate(p2, q)
    setInit!(body2; x = x2, q = q, F = F, τ = τ)
end

@inline getx3(body::Body, Δt) = getvnew(body) * Δt + body.x[2]
@inline getq3(body::Body, Δt) = Quaternion(Lmat(body.q[2]) * ωbar(body, Δt))
@inline getv1(body::Body, Δt) = (body.x[2] - body.x[1]) / Δt
@inline function getω1(body::Body, Δt)
    q1 = body.q[1]
    q2 = body.q[2]
    2 / Δt * (q1.s * imag(q2) - q2.s * imag(q1) - cross(imag(q1), imag(q2)))
    # 2/Δt*(VLᵀmat(body.q[1])*body.q[2])
end

# TODO use SOneTo and SUnitRange once they are faster
# @inline getvnew(body::Body) = body.s1[SOneTo(3)]
# @inline getωnew(body::Body) = body.s1[SUnitRange(4,6)]
@inline getvnew(body::Body) = body.s1[SVector(1, 2, 3)]
@inline getωnew(body::Body) = body.s1[SVector(4, 5, 6)]

@inline function derivωbar(body::Body{T}, Δt) where T
    ωnew = getωnew(body)
    msq = -sqrt(4 / Δt^2 - dot(ωnew, ωnew))
    Δt / 2 * [ωnew' / msq; SMatrix{3,3,T,9}(I)]
end

@inline function ωbar(body::Body, Δt)
    ωnew = getωnew(body)
    Δt / 2 * Quaternion(sqrt(4 / Δt^2 - dot(ωnew, ωnew)), ωnew)
end

@inline function dynamics(body::Body{T}, mechanism) where T
    No = mechanism.No
    Δt = mechanism.Δt

    ezg = SVector{3,T}(0, 0, -mechanism.g)
    dynT = body.m * ((getvnew(body) - getv1(body, Δt)) / Δt + ezg) - body.F[No]

    J = body.J
    ω1 = getω1(body, Δt)
    ωnew = getωnew(body)
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ωnew' * ωnew)
    dynR = skewplusdiag(ωnew, sq2) * (J * ωnew) - skewplusdiag(ω1, sq1) * (J * ω1) - 2 * body.τ[No]

    body.f = [dynT;dynR]

    for cid in connections(mechanism.graph, body.id)
        GtλTof!(body, geteqconstraint(mechanism, cid), mechanism)
    end

    for cid in ineqchildren(mechanism.graph, body.id)
        NtγTof!(body, getineqconstraint(mechanism, cid), mechanism)
    end

    return body.f
end

@inline function ∂dyn∂vel(body::Body{T}, Δt) where T
    J = body.J
    ωnew = getωnew(body)
    sq = sqrt(4 / Δt^2 - ωnew' * ωnew)

    dynT = SMatrix{3,3,T,9}(body.m / Δt * I)
    dynR = skewplusdiag(ωnew, sq) * J - J * ωnew * (ωnew' / sq) - skew(J * ωnew)

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end

@inline torqueFromForce(F::AbstractVector{T}, r::AbstractVector{T}) where T = cross(r,F)
@inline function setForce!(body::Body, F, τ, No)
    for i = 1:No
        body.F[i] = F
        body.τ[i] = τ
    end
end

