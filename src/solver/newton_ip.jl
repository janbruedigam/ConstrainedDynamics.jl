function newton_ip!(mechanism::Mechanism{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    ineqconstraints = mechanism.ineqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    for ineq in ineqconstraints
        ineq.s1 = ones(1)
        ineq.s0 = ones(1)
        ineq.γ1 = ones(1)
        ineq.γ0 = ones(1)
    end


    mechanism.μ = 1.0
    σ = 0.1

    meritf0 = meritf(mechanism)
    for n=Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph,ldu)
        solve!(graph,ldu,mechanism) # x̂1 for each body and constraint

        lineSearch!(mechanism,meritf0;iter=lineIter, warning=warning)


        meritf1 = meritf(mechanism)

        foreach(s1tos0!,bodies)
        foreach(s1tos0!,eqconstraints)
        foreach(s1tos0!,ineqconstraints)
        if normf(mechanism) < ε
            warning && display(n)
            return
        else
            while meritf1 < mechanism.μ
                mechanism.μ = σ*mechanism.μ
                meritf1 = meritf(mechanism)
            end
            meritf0 = meritf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf(mechanism),"."))
    end

    return
end

function lineSearch!(mechanism,meritf0;iter=10, warning::Bool=false)
    e = 0
    ldu = mechanism.ldu
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    ineqconstraints = mechanism.ineqconstraints

    computeα!(mechanism)
    αmax = mechanism.αmax

    for n=Base.OneTo(iter)
        for body in bodies
            lineStep!(body,getentry(ldu,body.id),e,αmax)# x1 = x0 - 1/(2^e)*d
        end
        for constraint in eqconstraints
            lineStep!(constraint,getentry(ldu,constraint.id),e,αmax)# x1 = x0 - 1/(2^e)*d
        end
        for constraint in ineqconstraints
            lineStep!(constraint,getineq(ldu,constraint.id),e,αmax)# x1 = x0 - 1/(2^e)*d
        end

        if meritf(mechanism) >= meritf0
            e += 1
        else
            return
        end
    end

    if warning
        display(string("WARNING:  lineSearch! did not converge. n = ",iter,"."))
    end
    return
end

@inline function lineStep!(node::Component,diagonal,e,α)
    node.s1 = node.s0 - 1/(2^e)*α*diagonal.Δs
    return
end

@inline function lineStep!(node::InequalityConstraint,entry,e,α)
    node.s1 = node.s0 - 1/(2^e)*α*entry.Δs
    node.γ1 = node.γ0 - 1/(2^e)*α*entry.Δγ
    return
end
