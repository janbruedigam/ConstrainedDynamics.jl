# xck = xdk+1
# vck = (xdk+1 - xdk)/Δt

# qck = qdk+1
# ωck = 2 V inv(qdk) (qdk+1 - qdk)/Δt
# ωckw = sqrt((2/Δt)^2 - ωckᵀωck) - 2/Δt

# Continuous values
@inline getx1(body::Body) = body.state.xc[1]
@inline getq1(body::Body) = body.state.qc[1]
@inline getx2(body::Body, Δt) = body.state.xc[1] + body.state.vc[2] * Δt
@inline getq2(body::Body, Δt) = body.state.qc[1] * ωbar(body.state.ωc[2], Δt)
@inline getv1(body::Body) = body.state.vc[1]
@inline getω1(body::Body) = body.state.ωc[1]
@inline getv2(body::Body) = body.state.vc[2]
@inline getω2(body::Body) = body.state.ωc[2]

# Discrete values
@inline getxd3(body::Body, Δt) = getx2(body, Δt)
@inline getqd3(body::Body, Δt) = getq2(body, Δt)
@inline getxd2(body::Body) = body.state.xd[2]
@inline getqd2(body::Body) = body.state.qd[2]
@inline getxd1(body::Body) = body.state.xd[1]
@inline getqd1(body::Body) = body.state.qd[1]

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
    xc = state.xc[1]
    qc = state.qc[1]
    vc = state.vc[1]
    ωc = state.ωc[1]

    state.xd[1] = xc - vc*Δt
    state.xd[2] = xc
    state.qd[1] = qc / ωbar(ωc,Δt)
    state.qd[2] = qc

    # just to set to some reasonable value
    state.xc[2] = xc + vc*Δt
    state.qc[2] = qc * ωbar(ωc,Δt)
    state.vc[2] = vc
    state.ωc[2] = ωc

    return
end

