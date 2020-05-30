@inline function settempvars!(body::Body{T}, x, v, F, q, ω, τ, f) where T
    state = body.state
    stateold = deepcopy(state)
    fold = body.f

    state.xc = x
    state.qc = q
    state.vc = v
    state.ωc = ω
    state.Fk[1] = F
    state.Fk[2] = F
    state.τk[1] = τ
    state.τk[2] = τ
    body.f = f
    return stateold, fold
end

@inline function resettempvars!(body::Body, state, f)
    body.state = state
    body.f = f
    return
end

@inline function ∂zp1∂z(mechanism, body::Body{T}, xc, vc, Fk, qc, ωc, τk, Δt) where T
    stateold, fold = settempvars!(body, xc, vc, Fk, qc, ωc, τk, zeros(T,6))


    state = body.state
    Z = @SMatrix zeros(T,3,3)
    Z2 = [Z Z]
    Z2t = Z2'
    E = SMatrix{3,3,T,9}(I)

    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    newton!(mechanism)

    # Position
    AposT = [E E*Δt]
    BposT = Z

        # This calculates the ϵ for q⊗Δq = q⊗(1 ϵᵀ)ᵀ
    AposR = VLmat(state.qk[1]) * [Rmat(ωbar(state.ωc, Δt))*LVᵀmat(state.qc) Lmat(state.qc)*derivωbar(state.ωc, Δt)]
    BposR = Z

    # Velocity
    AvelT = [Z E]
    BvelT = E*Δt/body.m

    J = body.J
    ω1 = state.ωc
    ω2 = state.ωsol[2]
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)

    ω1func = skewplusdiag(-ω1, sq1) * J - J * ω1 * (ω1' / sq1) + skew(J * ω1)
    ω2func = skewplusdiag(ω2, sq2) * J - J * ω2 * (ω2' / sq2) - skew(J * ω2)
    AvelR = [Z ω2func\ω1func]
    BvelR = 2*inv(ω2func)


    AT = [[AposT;AvelT] [Z2;Z2]]
    AR = [[Z2;Z2] [AposR;AvelR]]
    BT = [[BposT;BvelT] Z2t]
    BR = [Z2t [BposR;BvelR]]

    resettempvars!(body, stateold, fold)
    return [AT;AR], [BT;BR]
end