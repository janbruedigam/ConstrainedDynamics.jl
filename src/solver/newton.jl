# Newton with line search 
function newton!(mechanism::Mechanism{T,N,Nb,Ne,0}; ε = 1e-10, newtonIter = 100, lineIter = 10, warning::Bool = false) where {T,N,Nb,Ne}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    Δt = mechanism.Δt

    normf0 = normf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(mechanism) # x̂1 for each body and constraint

        normf1 = lineSearch!(mechanism, normf0;iter = lineIter, warning = warning)
        normΔs1 = normΔs(mechanism)

        foreachactive(updatesolution!, graph, bodies)
        foreachactive(updatesolution!, graph, eqcs)

        if normf1 < ε && normΔs1 < ε
            # warning && (@info string("Newton iterations: ",n))
            return
        else
            normf0 = normf1
        end
    end

    warning && (@info string("newton! did not converge. n = ", newtonIter, ", tol = ", normf(mechanism), "."))
    return
end

# Newton on interior point with line search
function newton!(mechanism::Mechanism{T,N,Nb,Ne,Ni}; ε = 1e-10, σ = 0.1, μ = 1.0, newtonIter = 100, lineIter = 10, warning::Bool = false) where {T,N,Nb,Ne,Ni}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    Δt = mechanism.Δt

    foreachactive(resetVars!, graph, ineqcs)
    mechanism.μ = μ

    meritf0 = meritf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(mechanism) # x̂1 for each body and constraint

        meritf1 = lineSearch!(mechanism, meritf0;iter = lineIter, warning = warning)

        foreachactive(updatesolution!, graph, bodies)
        foreachactive(updatesolution!, graph, eqcs)
        foreachactive(updatesolution!, graph, ineqcs)

        if normf(mechanism) < ε
            warning && (@info string("Newton iterations: ", n))
            return
        else
            while meritf1 < mechanism.μ
                mechanism.μ = σ * mechanism.μ
                meritf1 = meritf(mechanism)
            end
            meritf0 = meritf1
        end
    end

    warning && (@info string("newton! did not converge. n = ", newtonIter, ", tol = ", normf(mechanism), "."))
    return
end

function lineSearch!(mechanism::Mechanism{T,N,0}, normf0;iter = 10, warning::Bool = false) where {T,N}
    normf1 = normf0
    scale = 0
    ldu = mechanism.ldu
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for n = Base.OneTo(iter + 1)
        for body in bodies
            !body.active && continue
            lineStep!(body, getentry(ldu, body.id), scale)
        end
        for eqc in eqcs
            !eqc.active && continue
            lineStep!(eqc, getentry(ldu, eqc.id), scale)
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


function lineSearch!(mechanism::Mechanism, meritf0;iter = 10, warning::Bool = false)
    meritf1 = meritf0
    scale = 0
    ldu = mechanism.ldu
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints

    feasibilityStepLength!(mechanism)

    for n = Base.OneTo(iter)
        for body in bodies
            !body.active && continue
            lineStep!(body, getentry(ldu, body.id), scale, mechanism)
        end
        for eqc in eqcs
            !eqc.active && continue
            lineStep!(eqc, getentry(ldu, eqc.id), scale, mechanism)
        end
        for ineqc in ineqcs
            !ineqc.active && continue
            lineStep!(ineqc, getineqentry(ldu, ineqc.id), scale, mechanism)
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

@inline function lineStep!(body::Body, diagonal::DiagonalEntry, scale)
    body.state.vsol[2] = body.state.vsol[1] - 1 / (2^scale) * diagonal.Δs[SVector(1, 2, 3)]
    body.state.ωsol[2] = body.state.ωsol[1] - 1 / (2^scale) * diagonal.Δs[SVector(4, 5, 6)]
    return
end

@inline function lineStep!(eqc::EqualityConstraint, diagonal::DiagonalEntry, scale)
    eqc.λsol[2] = eqc.λsol[1] - 1 / (2^scale) * diagonal.Δs
    return
end

@inline function lineStep!(body::Body, diagonal, scale, mechanism)
    body.state.vsol[2] = body.state.vsol[1] - 1 / (2^scale) * mechanism.α * diagonal.Δs[SVector(1, 2, 3)]
    body.state.ωsol[2] = body.state.ωsol[1] - 1 / (2^scale) * mechanism.α * diagonal.Δs[SVector(4, 5, 6)]
    return
end

@inline function lineStep!(eqc::EqualityConstraint, diagonal, scale, mechanism)
    eqc.λsol[2] = eqc.λsol[1] - 1 / (2^scale) * mechanism.α * diagonal.Δs
    return
end

@inline function lineStep!(ineqc::InequalityConstraint, entry, scale, mechanism)
    ineqc.ssol[2] = ineqc.ssol[1] - 1 / (2^scale) * mechanism.α * entry.Δs
    ineqc.γsol[2] = ineqc.γsol[1] - 1 / (2^scale) * mechanism.α * entry.Δγ
    return
end