function newton_ip!(mechanism::Mechanism{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    ineqconstraints = mechanism.ineqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    for ineq in ineqconstraints
        ineq.sl1 = 1.
        ineq.sl0 = 1.
        ineq.ga1 = 1.
        ineq.ga0 = 1.
        ineq.slf1 = 1.
        ineq.slf0 = 1.
        ineq.psi1 = 1.
        ineq.psi0 = 1.
    end
    # for body in bodies
    #     body.s1 *= 1.
    #     body.s0 *= 1.
    # end
    # for constraint in eqconstraints
    #     constraint.s1 *= 0
    #     constraint.s0 *= 0
    # end

    # for ineq in ineqconstraints
    #     ineq.sl1 *= 100
    #     ineq.sl0 *= 100
    #     ineq.ga1 *= 100
    #     ineq.ga0 *= 100
    # end
    #
    # slgan = 0.0
    # for ineq in mechanism.ineqconstraints
    #     slgan += ineq.sl1*ineq.ga1
    # end

    mechanism.μ = 1.#slgan/length(ineqconstraints)
    σ = 0.1

    # normf0 = normf(mechanism)
    meritf0 = meritf(mechanism)
    for n=Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph,ldu)
        solve!(graph,ldu,mechanism) # x̂1 for each body and constraint

        lineSearch!(mechanism,meritf0;iter=lineIter, warning=warning)

        # foreach(update!,bodies,ldu,αsmax)
        # foreach(update!,eqconstraints,ldu,αγmax)
        # foreach(update!,ineqconstraints,ldu,αsmax,αγmax)

        meritf1 = meritf(mechanism)

        normsol = normΔs(mechanism)
        foreach(s1tos0!,bodies)
        foreach(s1tos0!,eqconstraints)
        foreach(s1tos0!,ineqconstraints)
        if normf(mechanism) < ε # && normsol < ε
            display(n)
            return
        else
            if meritf1 < mechanism.μ*0.01 #&& mechanism.μ > ε*0.1
                mechanism.μ = σ*mechanism.μ
                meritf0 = meritf(mechanism)
            else
                meritf0=meritf1
            end
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
    node.s1 = node.s0 - 1/(2^e)*α*diagonal.ŝ
    return
end

@inline function lineStep!(node::InequalityConstraint,entry,e,α)
    node.sl1 = node.sl0 - 1/(2^e)*α*entry.sl
    node.ga1 = node.ga0 - 1/(2^e)*α*entry.ga
    node.slf1 = node.slf0 - 1/(2^e)*α*entry.slf
    node.psi1 = node.psi0 - 1/(2^e)*α*entry.psi
    node.b1 = node.b0 - 1/(2^e)*α*entry.b
    return
end
