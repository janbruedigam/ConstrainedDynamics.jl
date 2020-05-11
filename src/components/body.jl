mutable struct Body{T} <: AbstractBody{T}
    id::Int64
    name::String

    m::T
    J::SMatrix{3,3,T,9}

    x::Vector{SVector{3,T}}
    q::Vector{Quaternion{T}}

    F::Vector{SVector{3,T}}
    τ::Vector{SVector{3,T}}

    s0::SVector{6,T}
    s1::SVector{6,T}
    f::SVector{6,T}


    function Body(m::T, J::AbstractArray{T,2}; name::String="") where T
        x = [zeros(T, 3)]
        q = [Quaternion{T}()]

        F = [zeros(T, 3)]
        τ = [zeros(T, 3)]

        s0 = zeros(T, 6)
        s1 = zeros(T, 6)
        f = zeros(T, 6)

        new{T}(getGlobalID(), name, m, J, x, q, F, τ, s0, s1, f)
    end

    function Body(shape::Shape; name::String="")
        body = Body(shape.m, shape.J; name=name)
        push!(shape.bodyids, body.id)
        return body
    end

    function Body(::Type{T},contents...) where T
        new{T}(getGlobalID(), contents...)
    end
end

mutable struct Origin{T} <: AbstractBody{T}
    id::Int64
    name::String

    Origin{T}(; name::String="") where T = new{T}(getGlobalID(), name)

    function Origin(shape::Shape{T}; name::String="") where T
        origin = Origin{T}(name=name)
        push!(shape.bodyids, origin.id)
        return origin
    end
end

function Base.deepcopy(b::Body{T}) where T
    contents = []
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    Body(T, contents...)
end

function Base.deepcopy(b::Body{T}, shape::Shape) where T
    contents = []
    
    for i = 2:getfieldnumber(b)
        push!(contents, deepcopy(getfield(b, i)))
    end

    body = Body(T, contents...)
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

@inline getM(body::Body{T}) where T = [[SMatrix{3,3,T,9}(I*body.m);@SMatrix zeros(T,3,3)] [@SMatrix zeros(T,3,3);body.J]]

@inline getx3(body::Body, Δt) = getvnew(body) * Δt + body.x[2]
@inline getq3(body::Body, Δt) = Quaternion(Lmat(body.q[2]) * ωbar(body, Δt))
#TODO why not simple from s0?
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

@inline function settempvars(body::Body{T}, x, v, F, q, ω, τ, s0, s1, f, Δt) where T
    xold = body.x[2]
    vold = getv1(body,Δt)
    Fold = body.F[2]
    qold = body.q[2]
    ωold = getω1(body,Δt)
    τold = body.τ[2]
    s0old = body.s0
    s1old = body.s1
    fold = body.f
    

    body.x[2] = x
    body.x[1] = x - v*Δt
    body.F[2] = F
    body.q[2] = q
    body.q[1] = Δt/2 * (q * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), -SVector{3,T}(ω))) # this accounts for the non-unit-quaternion inverse (q2 * conj(ωbar))
    body.τ[2] = τ
    body.s0 = s0
    body.s1 = s1
    body.f = f
    return xold, vold, Fold, qold, ωold, τold, s0old, s1old, fold
end

@inline function dynamics(mechanism, body::Body{T}) where T
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
        GtλTof!(mechanism, body, geteqconstraint(mechanism, cid))
    end

    for cid in ineqchildren(mechanism.graph, body.id)
        NtγTof!(mechanism, body, getineqconstraint(mechanism, cid))
    end

    return body.f
end

@inline function ∂dyn∂pos(body::Body{T}, Δt) where T
    J = body.J
    ωnew = getωnew(body)
    sq = sqrt(4 / Δt^2 - ωnew' * ωnew)

    dynT = SMatrix{3,3,T,9}(body.m / Δt^2 * I)
    dynR = (skewplusdiag(ωnew, sq) * J - J * ωnew * (ωnew' / sq) - skew(J * ωnew)) * 2/Δt * VLᵀmat(body.q[2])*LVᵀmat(getq3(body,Δt))

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
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

# TODO tensor terms
@inline function ∂zp1∂z(mechanism, body::Body{T}, xd, vd, Fd, qd, ωd, τd, Δt) where T
    xold, vold, Fold, qold, ωold, τold, s0old, s1old, fold = settempvars(body, xd, vd, Fd, qd, ωd, τd, [vd;ωd], [vd;ωd], zeros(6,T), Δt)

    newton!(mechanism)

    # Velocity
    Z = @SMatrix zeros(T,3,3)
    Z43 = @SMatrix zeros(T,4,3)
    Z2 = [Z Z]
    Z432 = [Z43 Z43]
    E = SMatrix{3,3,T,9}(I)

    AvelT = [Z E]
    BvelT = E*Δt/body.m

    J = body.J
    ω1 = getω1(body, Δt)
    ωnew = getωnew(body)
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ωnew' * ωnew)
    dynR = skewplusdiag(ωnew, sq2) * (J * ωnew) - skewplusdiag(ω1, sq1) * (J * ω1) - 2 * body.τ[No]

    ω1func = skewplusdiag(-ω1, sq1) * J - J * ω1 * (ω1' / sq1) + skew(J * ω1)
    ω2func = skewplusdiag(ωnew, sq2) * J - J * ωnew * (ωnew' / sq2) - skew(J * ωnew)
    AvelR = [Z ω2func\ω1func]
    BvelR = 2*E

    # Position
    AposT = [E Z] + AvelT*Δt
    BposT = BvelT*Δt

    AposR = Δt/2*([Lmat(ωbar(body, Δt))*LVᵀmat(qd) Z] + Lmat(qd)*derivωbar(body, Δt)*AvelR)
    BposR = Δt/2*Lmat(qd)*derivωbar(body, Δt)*BvelR

    AT = [[AposT;AvelT] [Z2;Z2]]
    AR = [[Z432;Z2] [AposR;AvelR]]
    BT = [BposT;BvelT]
    BR = [BposR;BvelR]

    settempvars(body, xold, vold, Fold, qold, ωold, τold, s0old, s1old, fold, Δt)
    return [AT;AR], [BT;BR]
end


# @inline function ∂dyn∂posm1(body::Body{T}, Δt) where T
#     J = body.J
#     ω1 = getω1(body, Δt)
#     sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
#     ωnew = getωnew(body)
#     sqnew = sqrt(4 / Δt^2 - ωnew' * ωnew)

#     dynT = SMatrix{3,3,T,9}(-2*body.m / Δt^2 * I)
#     dynR = -(skewplusdiag(ωnew, sqnew) * J - J * ωnew * (ωnew' / sqnew) - skew(J * ωnew)) * 2/Δt * VRmat(getq3(body,Δt))*RᵀVᵀmat(body.q[2]) +
#         (skewplusdiag(ω1, -sq1) * J + J * ω1 * (ω1' / sq1) - skew(J * ω1)) * 2/Δt * VLᵀmat(body.q[1])*LVᵀmat(body.q[2])

#     Z = @SMatrix zeros(T, 3, 3)

#     return [[dynT; Z] [Z; dynR]]
# end

# @inline function ∂dyn∂velm1(body::Body{T}, Δt) where T
#     J = body.J
#     ω1 = getω1(body, Δt)
#     sq = sqrt(4 / Δt^2 - ω1' * ω1)

#     dynT = SMatrix{3,3,T,9}(-body.m / Δt * I)
#     dynR = skewplusdiag(ω1, -sq) * J + J * ω1 * (ω1' / sq) - skew(J * ω1)

#     Z = @SMatrix zeros(T, 3, 3)

#     return [[dynT; Z] [Z; dynR]]
# end

# @inline function ∂dyn∂con(body::Body{T}, Δt) where T
#     dynT = SMatrix{3,3,T,9}(-I)
#     dynR = SMatrix{3,3,T,9}(-2*I)

#     Z = @SMatrix zeros(T, 3, 3)

#     return [[dynT; Z] [Z; dynR]]
# end

@inline torqueFromForce(F::AbstractVector{T}, r::AbstractVector{T}) where T = cross(r, F)
@inline function setForce!(body::Body, F, τ, No)
    for i = 1:No
        body.F[i] = F
        body.τ[i] = τ
    end
end

