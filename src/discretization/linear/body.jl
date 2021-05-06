@inline function dynamics(mechanism, body::Body{T}) where T
    state = body.state
    Δt = mechanism.Δt
    graph = mechanism.graph

    ezg = SA{T}[0; 0; -mechanism.g]
    dynT = body.m * ((state.vsol[2] - state.vc) / Δt + ezg) - state.Fk[1]

    J = body.J
    ω1 = state.ωc
    ω2 = state.ωsol[2]
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    dynR = skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1) - 2 * state.τk[1]

    state.d = [dynT;dynR]

    for connectionid in connections(graph, body.id)
        GtλTof!(mechanism, body, geteqconstraint(mechanism, connectionid))
    end

    for connectionid in springconnections(graph, body.id)
        # display("here")
        springTof!(mechanism, body, geteqconstraint(mechanism, connectionid))
    end

    for connectionid in damperconnections(graph, body.id)
        damperTof!(mechanism, body, geteqconstraint(mechanism, connectionid))
    end

    for childid in ineqchildren(graph, body.id)
        NtγTof!(mechanism, body, getineqconstraint(mechanism, childid))
    end

    return state.d
end
