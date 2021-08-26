# Newton with line search 
function newton!(mechanism::Mechanism{T,Nn,Ne,Nb,Nf,0};
        ε = 1e-10, newtonIter = 100, lineIter = 10, warning::Bool = false
    ) where {T,Nn,Ne,Nb,Nf}
    
    system = mechanism.system
    eqcs = mechanism.eqconstraints
    bodies = mechanism.bodies

    normf0 = normf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        ldu_solve!(system) # x̂1 for each body and constraint

        normf1 = lineSearch!(mechanism, normf0;iter = lineIter, warning = warning)
        normΔs1 = normΔs(mechanism)

        foreach(updatesolution!, bodies)
        foreach(updatesolution!, eqcs)

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
function newton!(mechanism::Mechanism;
        ε = 1e-10, σ = 0.1, μ = 1.0, newtonIter = 100, lineIter = 10, warning::Bool = false
    )
    
    system = mechanism.system
    eqcs = mechanism.eqconstraints
    bodies = mechanism.bodies
    frics = mechanism.frictions
    ineqcs = mechanism.ineqconstraints

    foreach(resetVars!, ineqcs)
    mechanism.μ = μ

    meritf0 = meritf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        ldu_solve!(system) # x̂1 for each body and constraint

        meritf1 = lineSearch!(mechanism, meritf0;iter = lineIter, warning = warning)

        foreach(updatesolution!, bodies)
        foreach(updatesolution!, eqcs)
        foreach(updatesolution!, ineqcs)
        foreach(updatesolution!, frics)

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
