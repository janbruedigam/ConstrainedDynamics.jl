# xck = xdk+1
# vck = (xdk+1 - xdk)/Δt

# qck = qdk+1
# ωck = 2 V inv(qdk) (qdk+1 - qdk)/Δt
# ωckw = sqrt((2/Δt)^2 - ωckᵀωck) - 2/Δt

# Continuous values
@inline getx1(body::Body) = body.state.xd[2]
@inline getq1(body::Body) = body.state.qd[2]
@inline getx2(body::Body, Δt) = body.state.xd[2] + getv2(body) * Δt
@inline getq2(body::Body, Δt) = Quaternion(Lmat(body.state.qd[2]) * ωbar(body.state.ωc[2], Δt))
@inline getv1(body::Body) = body.state.vc[1]
@inline getω1(body::Body) = body.state.ωc[1]
@inline getv2(body::Body) = body.state.vc[2]
@inline getω2(body::Body) = body.state.ωc[2]

@inline function derivωbar(ω::SVector{3,T}, Δt) where T
    msq = -sqrt(4 / Δt^2 - dot(ω, ω))
    Δt / 2 * [ω' / msq; SMatrix{3,3,T,9}(I)]
end

@inline function ωbar(ω, Δt)
    Δt / 2 * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), ω)
end

@inline function setForce!(body::Body, F, τ)
    body.F[2] = F
    body.τ[2] = τ
end


@inline function discretizestate!(body::Body, Δt)
    state = body.state
    xd = state.xc[1]
    qd = state.qc[1]
    vc = state.vc[1]
    ωc = state.ωc[1]

    state.xd[1] = xd - vc*Δt
    state.xd[2] = xd
    state.qd[1] = qd / ωbar(ωc,Δt)
    state.qd[2] = qd

    # just to set to some reasonable value
    state.vc[2] = vc
    state.ωc[2] = ωc

    return
end


@inline function dynamics(mechanism, body::Body{T}) where T
    Δt = mechanism.Δt

    ezg = SVector{3,T}(0, 0, -mechanism.g)
    dynT = body.m * ((getv2(body) - getv1(body)) / Δt + ezg) - body.F[2]

    J = body.J
    ω1 = getω1(body)
    ω2 = getω2(body)
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    dynR = skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1) - 2 * body.τ[2]

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
    ω2 = getω2(body)
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    dynT = SMatrix{3,3,T,9}(body.m / Δt^2 * I)
    dynR = (skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2)) * 2/Δt * VLᵀmat(body.state.qd[2])*LVᵀmat(getq2(body,Δt))

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end

@inline function ∂dyn∂vel(body::Body{T}, Δt) where T
    J = body.J
    ω2 = getω2(body)
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    dynT = SMatrix{3,3,T,9}(body.m / Δt * I)
    dynR = skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2)

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end