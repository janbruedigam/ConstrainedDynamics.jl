@inline function dynamics(mechanism, body::Body{T}) where T
    state = body.state
    Δt = mechanism.Δt
    graph = mechanism.graph

    ezg = SVector{3,T}(0, 0, -mechanism.g)
    dynT = body.m * ((state.vsol[2] - state.vc) / Δt + ezg) - state.Fk[1]

    J = body.J
    ω1 = state.ωc
    ω2 = state.ωsol[2]
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    dynR = skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1) - 2 * state.τk[1]

    state.d = [dynT;dynR]

    for childid in connections(mechanism.graph, body.id)
        isinactive(graph, childid) && continue
        GtλTof!(mechanism, body, geteqconstraint(mechanism, childid))
    end

    for childid in ineqchildren(mechanism.graph, body.id)
        isinactive(graph, childid) && continue
        NtγTof!(mechanism, body, getineqconstraint(mechanism, childid))
    end

    return state.d
end

@inline function ∂dyn∂vel(body::Body{T}, Δt) where T
    J = body.J
    ω2 = body.state.ωsol[2]
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    dynT = SMatrix{3,3,T,9}(body.m / Δt * I)
    dynR = skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2)

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end

@inline function ∂zp1∂z(mechanism, body::Body{T}, Δt) where T
    state = body.state
    Z = @SMatrix zeros(T,3,3)
    Z2 = [Z Z]
    Z2t = Z2'
    E = SMatrix{3,3,T,9}(I)

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

    return [AT;AR], [BT;BR]
end

