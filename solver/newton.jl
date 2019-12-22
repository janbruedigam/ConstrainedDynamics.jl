function newton!(robot::Robot{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    # n = 1
    links = robot.links
    constraints = robot.constraints
    graph = robot.graph
    ldu = robot.ldu

    normf0 = normf(robot)
    for n=Base.OneTo(newtonIter)
        setentries!(robot)
        factor!(graph,ldu)
        solve!(graph,ldu) # x̂1 for each link and constraint
        correctλ!(robot)

        foreach(update!,links,ldu)
        foreach(update!,constraints,ldu)

        normf1 = normf(robot)
        normf1>normf0 && lineSearch!(robot,normf0;iter=lineIter, warning=warning)

        if normΔs(robot) < ε && normf1 < ε
            foreach(s1tos0!,links)
            foreach(s1tos0!,constraints)
            # display(n)
            return
        else
            foreach(s1tos0!,links)
            foreach(s1tos0!,constraints)
            normf0=normf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf0,"."))
    end

    return
end

function lineSearch!(robot,normf0;iter=10, warning::Bool=false)
    α = 1
    ldu = robot.ldu
    links = robot.links
    constraints = robot.constraints
    for link in links
        lineStep!(link,getentry(ldu,link.id),α)# x1 = x0 + 1/(2^α)*d
    end
    for constraint in constraints
        lineStep!(constraint,getentry(ldu,constraint.id),α)# x1 = x0 + 1/(2^α)*d
    end

    for n=1:iter
        α += 1
        if normf(robot) >= normf0
            for link in links
                lineStep!(link,getentry(ldu,link.id),α)# x1 = x0 + 1/(2^α)*d
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
