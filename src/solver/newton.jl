function newton!(mechanism::Mechanism{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    # n = 1
    bodies = mechanism.bodies
    constraints = mechanism.eqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    normf0 = normf(mechanism)
    for n=Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph,ldu)
        solve!(graph,ldu,mechanism) # x̂1 for each body and constraint

        foreach(update!,bodies,ldu)
        foreach(update!,constraints,ldu)
        # foreach(checkωnorm!,bodies,dt)

        normf1 = normf(mechanism)
        normf1>normf0 && lineSearch!(mechanism,normf0;iter=lineIter, warning=warning)

        # Move foreach out of if-else? Only if calculating normΔs earlier !!!
        if normΔs(mechanism) < ε && normf1 < ε
            foreach(s1tos0!,bodies)
            foreach(s1tos0!,constraints)
            # display(n)
            return
        else
            foreach(s1tos0!,bodies)
            foreach(s1tos0!,constraints)
            normf0=normf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf0,"."))
    end

    return
end

function lineSearch!(mechanism,normf0;iter=10, warning::Bool=false)
    α = 1
    ldu = mechanism.ldu
    bodies = mechanism.bodies
    constraints = mechanism.constraints
    for body in bodies
        lineStep!(body,getentry(ldu,body.id),α)# x1 = x0 + 1/(2^α)*d
    end
    for constraint in constraints
        lineStep!(constraint,getentry(ldu,constraint.id),α)# x1 = x0 + 1/(2^α)*d
    end

    for n=Base.OneTo(iter)
        α += 1
        if normf(mechanism) >= normf0
            for body in bodies
                lineStep!(body,getentry(ldu,body.id),α)# x1 = x0 + 1/(2^α)*d
            end
            for constraint in constraints
                lineStep!(constraint,getentry(ldu,constraint.id),α)# x1 = x0 + 1/(2^α)*d
            end
        else
            return
        end
    end

    if warning
        display(string("WARNING:  lineSearch! did not converge. n = ",iter,"."))
    end
    return
end

@inline lineStep!(node,diagonal,α) = (node.s1 = node.s0 - 1/(2^α)*diagonal.ŝ; return)
