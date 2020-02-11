function newton_ip!(mechanism::Mechanism{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    # n = 1
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    ineqconstraints = mechanism.ineqconstraints
    graph = mechanism.graph
    ldu = mechanism.ldu
    dt = mechanism.dt

    for ineq in mechanism.ineqconstraints
        ineq.sl1 = rand()
        ineq.sl0 = rand()
        ineq.ga1 = rand()
        ineq.ga0 = rand()
    end

    mechanism.μ = 1.
    σ = 0.2

    normf0 = normf(mechanism)
    for n=Base.OneTo(newtonIter)
        setentries!(mechanism)
        factor!(graph,ldu)
        solve!(graph,ldu,mechanism) # x̂1 for each body and constraint

        # for ineqs in mechanism.ineqconstraints
        #     body = getbody(mechanism,ineqs.pid)
        #     computeΔ!(body,getentry(ldu,body.id),mechanism)
        # end
        computeα!(mechanism)

        αsmax = mechanism.αsmax
        αγmax = mechanism.αγmax

        foreach(update!,bodies,ldu,αsmax)
        foreach(update!,eqconstraints,ldu,αγmax)
        foreach(update!,ineqconstraints,ldu,αsmax,αγmax)

        normf1 = normf(mechanism)
        # normf1>normf0 && lineSearch!(mechanism,normf0;iter=lineIter, warning=warning)

        # normΔs not changed yet !!!
        if normΔs(mechanism) < ε && normf1 < ε
            foreach(s1tos0!,bodies)
            foreach(s1tos0!,eqconstraints)
            foreach(s1tos0!,ineqconstraints)
            # display(n)
            return
        else
            foreach(s1tos0!,bodies)
            foreach(s1tos0!,eqconstraints)
            foreach(s1tos0!,ineqconstraints)
            mechanism.μ = σ*mechanism.μ
            normf0=normf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf0,"."))
    end

    return
end
