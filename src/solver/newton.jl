# Newton with line search 
function newton!(mechanism::Mechanism{T,N,Nb,Ne,0};
        ε = 1e-10, newtonIter = 100, lineIter = 10, warning::Bool = false
    ) where {T,N,Nb,Ne}
    
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    normf0 = normf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(mechanism) # x̂1 for each body and constraint

        normf1 = lineSearch!(mechanism, normf0;iter = lineIter, warning = warning)
        normΔs1 = normΔs(mechanism)

        foreachactive(updatesolution!, bodies)
        foreachactive(updatesolution!, eqcs)

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
function newton!(mechanism::Mechanism{T,N,Nb,Ne,Ni};
        ε = 1e-10, σ = 0.1, μ = 1.0, newtonIter = 100, lineIter = 10, warning::Bool = false
    ) where {T,N,Nb,Ne,Ni}
    
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu

    foreachactive(resetVars!, ineqcs)
    mechanism.μ = μ

    meritf0 = meritf(mechanism)
    for n = Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph, ldu)
        solve!(mechanism) # x̂1 for each body and constraint

        meritf1 = lineSearch!(mechanism, meritf0;iter = lineIter, warning = warning)

        foreachactive(updatesolution!, bodies)
        foreachactive(updatesolution!, eqcs)
        foreachactive(updatesolution!, ineqcs)

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

function newton!(mechanism::LinearMechanism{T,N,Nb,Ne,0};
        ε = 1e-10, newtonIter = 100, lineIter = 10, warning::Bool = false
    ) where {T,N,Nb,Ne}

    A = mechanism.A
    Bu = mechanism.Bu
    Bλ = mechanism.Bλ
    G = mechanism.G
    nc = size(G)[1]
    Z = zeros(T,nc,nc)

    normf0 = normf(mechanism)
    for n = Base.OneTo(newtonIter)
        M = [-I Bλ;G Z]
        Δsol = M\[fdynamics(mechanism); fconstraints(mechanism)]
        mechanism.Δz = Δsol[1:12*Nb]
        mechanism.Δλ = Δsol[12*Nb+1:end]

        normf1 = lineSearch!(mechanism, normf0;iter = lineIter, warning = warning)
        # normΔs1 = normΔs(mechanism)

        mechanism.zsol[1] = mechanism.zsol[2]
        mechanism.λsol[1] = mechanism.λsol[2]

        if normf1 < ε #&& normΔs1 < ε
            # warning && (@info string("Newton iterations: ",n))
            return
        else
            normf0 = normf1
        end
    end

    warning && (@info string("newton! did not converge. n = ", newtonIter, ", tol = ", normf(mechanism), "."))
    return
end
