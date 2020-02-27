# Newton with line search 
function newton!(mechanism::Mechanism{T,Nl,0}; ε = 1e-10, newtonIter = 100, lineIter = 10, warning::Bool = false) where {T,Nl}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    normf0 = normf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(graph, ldu, mechanism) # x̂1 for each body and constraint

        normf1 = lineSearch!(mechanism, normf0;iter = lineIter, warning = warning)
        normΔs1 = normΔs(mechanism)

        foreach(s1tos0!, bodies)
        foreach(s1tos0!, eqcs)

        if normΔs1 < ε && normf1 < ε
            # warning && (@info string("Newton iterations: ",n))
            return
        else
            normf0 = normf1
        end
    end

    warning && (@info string("newton_ip! did not converge. n = ", newtonIter, ", tol = ", normf(mechanism), "."))
    return
end

# Newton on interior point with line search
function newton!(mechanism::Mechanism{T,Nl}; ε = 1e-10, σ = 0.1, newtonIter = 100, lineIter = 10, warning::Bool = false) where {T,Nl}
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    foreach(resetVars!, ineqcs)

    meritf0 = meritf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(graph, ldu, mechanism) # x̂1 for each body and constraint

        meritf1 = lineSearch!(mechanism, meritf0;iter = lineIter, warning = warning)

        foreach(s1tos0!, bodies)
        foreach(s1tos0!, eqcs)
        foreach(s1tos0!, ineqcs)

        if normf(mechanism) < ε
            warning && (@info string("Newton iterations: ", n))
            return
        else
            for ineq in ineqcs
                if meritf1 < ineq.μ
                    ineq.μ = σ * ineq.μ
                    meritf1 = meritf(mechanism)
                end
                # while meritf1 < ineq.μ
                #     ineq.μ = σ*ineq.μ
                #     meritf1 = meritf(mechanism)
                # end
            end
            meritf0 = meritf1
        end
    end

    warning && (@info string("newton_ip! did not converge. n = ", newtonIter, ", tol = ", normf(mechanism), "."))
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
            lineStep!(body, getentry(ldu, body.id), scale)
        end
        for eqc in eqcs
            lineStep!(eqc, getentry(ldu, eqc.id), scale)
        end

        normf1 = normf(mechanism)
        if normf1 >= normf0
            scale += 1
        else
            return normf1
        end
    end

    warning && (@info string("lineSearch! did not converge. n = ", iter, "."))
    return normf1
end


function lineSearch!(mechanism::Mechanism, meritf0;iter = 10, warning::Bool = false)
    meritf1 = meritf0
    scale = 0
    ldu = mechanism.ldu
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints

    α = 1.
    for ineqc in ineqcs
        αmax = feasibilityStepLenth(ineqc, ldu)
        (αmax < α) && (α = αmax)
    end

    for n = Base.OneTo(iter)
        for body in bodies
            lineStep!(body, getentry(ldu, body.id), scale, α)
        end
        for eqc in eqcs
            lineStep!(eqc, getentry(ldu, eqc.id), scale, α)
        end
        for ineqc in ineqcs
            lineStep!(ineqc, getineqentry(ldu, ineqc.id), scale, α)
        end

        meritf1 = meritf(mechanism)
        if meritf1 >= meritf0
            scale += 1
        else
            return meritf1
        end
    end

    warning && (@info string("lineSearch! did not converge. n = ", iter, "."))
    return meritf1
end

@inline function lineStep!(component::Component, diagonal::DiagonalEntry, scale)
    component.s1 = component.s0 - 1 / (2^scale) * diagonal.Δs
    return
end

@inline function lineStep!(node::Component, diagonal, scale, α)
    node.s1 = node.s0 - 1 / (2^scale) * α * diagonal.Δs
    return
end

@inline function lineStep!(ineqc::InequalityConstraint, entry, scale, α)
    ineqc.s1 = ineqc.s0 - 1 / (2^scale) * α * entry.Δs
    ineqc.γ1 = ineqc.γ0 - 1 / (2^scale) * α * entry.Δγ
    return
end