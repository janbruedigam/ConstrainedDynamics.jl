# L(xck,vck) -> 1/2 Δt (Ld(xdk,(xdk+1-xdk)/Δt) + Ld(xdk+1,(xdk+1-xdk)/Δt))
# L(qck,ωck) -> 1/2 Δt (Ld(qdk,2 V qdk† (qdk+1-qdk)/Δt) + Ld(qdk+1,2 V qdk† (qdk+1-qdk)/Δt)
# ωckw = sqrt((2/Δt)^2 - ωckᵀωck) - 2/Δt
# 1/2 (Fdk + Fdk+1)

METHODORDER = 2
getGlobalOrder() = (global METHODORDER; return METHODORDER)

# Convenience functions
@inline getx3(state::State, Δt) = state.xk[2] + state.vsol[2]*Δt
@inline getq3(state::State, Δt) = state.qk[2] * ωbar(state.ωsol[2],Δt)

@inline function derivωbar(ω::SVector{3,T}, Δt) where T
    msq = -sqrt(4 / Δt^2 - dot(ω, ω))
    Δt / 2 * [ω' / msq; SMatrix{3,3,T,9}(I)]
end

@inline function ωbar(ω, Δt)
    Δt / 2 * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), ω)
end

@inline function setForce!(state::State, F, τ)
    state.Fk[2] = F
    state.τk[2] = τ
    return
end


@inline function discretizestate!(body::Body, Δt)
    state = body.state
    xc = state.xc
    qc = state.qc
    vc = state.vc
    ωc = state.ωc

    state.xk[1] = xc
    state.xk[2] = xc + vc*Δt
    state.qk[1] = qc
    state.qk[2] = qc * ωbar(ωc,Δt)

    return
end

@inline function updatestate!(body::Body{T}, Δt) where T
    state = body.state

    state.xc = state.xsol[2]
    state.qc = state.qsol[2]
    state.vc = state.vsol[2]
    state.ωc = state.ωsol[2]

    state.xk[1] = state.xk[2]
    state.xk[2] = state.xk[2] + state.vsol[2]*Δt
    state.qk[1] = state.qk[2]
    state.qk[2] = state.qk[2] * ωbar(state.ωsol[2],Δt)

    state.xsol[2] = state.xk[2]
    state.qsol[2] = state.qk[2]

    state.Fk[1] = state.Fk[2]
    state.Fk[2] = @SVector zeros(T,3)
    state.τk[1] = state.τk[2]
    state.τk[2] = @SVector zeros(T,3)
    return
end

@inline function setsolution!(body::Body)
    state = body.state
    state.xsol[2] = state.xk[2]
    state.qsol[2] = state.qk[2]
    state.vsol[1] = state.vc
    state.vsol[2] = state.vc
    state.ωsol[1] = state.ωc
    state.ωsol[2] = state.ωc
    return
end