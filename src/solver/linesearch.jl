function lineSearch!(mechanism::Mechanism{T,Nn,Nb,Ne,0}, normf0; iter = 10, warning::Bool = false) where {T,Nn,Nb,Ne}
    normf1 = normf0
    scale = 0
    system = mechanism.system
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for n = Base.OneTo(iter + 1)
        for body in bodies
            lineStep!(body, getentry(system, body.id), scale)
            if norm(body.state.ωsol[2]) > 1/mechanism.Δt
                error("Excessive angular velocity. Body-ID: "*string(body.id)*", ω: "*string(body.state.ωsol[2])*".")
            end
        end
        for eqc in eqcs
            lineStep!(eqc, getentry(system, eqc.id), scale)
        end

        normf1 = normf(mechanism)
        if normf1 >= normf0
            scale += 1
        else
            return normf1
        end
    end

    warning && (@info string("lineSearch! did not converge. n = ", iter, ". Last tol: ", normf1))
    return normf1
end

function lineSearch!(mechanism::Mechanism, meritf0; iter = 10, warning::Bool = false)
    meritf1 = meritf0
    scale = 0
    system = mechanism.system
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints

    feasibilityStepLength!(mechanism)

    for n = Base.OneTo(iter)
        for body in bodies
            lineStep!(body, getentry(system, body.id), scale, mechanism)
            if norm(body.state.ωsol[2]) > 1/mechanism.Δt
                error("Excessive angular velocity. Body-ID: "*string(body.id)*", ω: "*string(body.state.ωsol[2])*".")
            end
        end
        for eqc in eqcs
            lineStep!(eqc, getentry(system, eqc.id), scale, mechanism)
        end
        for ineqc in ineqcs
            lineStep!(ineqc, getentry(system, ineqc.id), scale, mechanism)
        end

        meritf1 = meritf(mechanism)
        if meritf1 >= meritf0
            scale += 1
        else
            return meritf1
        end
    end

    warning && (@info string("lineSearch! did not converge. n = ", iter, ". Last tol: ", meritf1))
    return meritf1
end


@inline function lineStep!(body::Body, vector_entry::Entry, scale)
    body.state.vsol[2] = body.state.vsol[1] + 1 / (2^scale) * vector_entry.value[SA[1; 2; 3]]
    body.state.ωsol[2] = body.state.ωsol[1] + 1 / (2^scale) * vector_entry.value[SA[4; 5; 6]]
    return
end

@inline function lineStep!(eqc::EqualityConstraint, vector_entry::Entry, scale)
    eqc.λsol[2] = eqc.λsol[1] + 1 / (2^scale) * vector_entry.value
    return
end

@inline function lineStep!(body::Body, vector_entry::Entry, scale, mechanism)
    body.state.vsol[2] = body.state.vsol[1] + 1 / (2^scale) * mechanism.α * vector_entry.value[SA[1; 2; 3]]
    body.state.ωsol[2] = body.state.ωsol[1] + 1 / (2^scale) * mechanism.α * vector_entry.value[SA[4; 5; 6]]
    return
end

@inline function lineStep!(eqc::EqualityConstraint, vector_entry::Entry, scale, mechanism)
    eqc.λsol[2] = eqc.λsol[1] + 1 / (2^scale) * mechanism.α * vector_entry.value
    return
end

@inline function lineStep!(ineqc::InequalityConstraint{T,N}, vector_entry::Entry, scale, mechanism) where {T,N}
    ineqc.ssol[2] = ineqc.ssol[1] + 1 / (2^scale) * mechanism.α * vector_entry.value[SVector{N,Int64}(1:N)]
    ineqc.γsol[2] = ineqc.γsol[1] + 1 / (2^scale) * mechanism.α * vector_entry.value[SVector{N,Int64}(N+1:2*N)]
    return
end
