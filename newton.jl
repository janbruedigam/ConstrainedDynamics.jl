function newton!(robot::Robot{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=20, warning::Bool=false) where {T,Nl}
    n = 1
    links = robot.links
    constraints = robot.constraints
    graph = robot.graph
    ldu = robot.ldu

    normf0 = normf(robot)
    for outer n=1:newtonIter
        setentries!(robot)
        factor!(graph,ldu)
        solve!(graph,ldu) # x̂1 for each link and constraint
        # correctλ!(robot) # for simplified system
        for link in links
            update!(link,getentry(ldu,link.id))
        end
        for constraint in constraints
            update!(constraint,getentry(ldu,constraint.id))
        end

        normf1 = normf(robot)
        normf1>normf0 ? lineSearch!(robot,normf0;iter=lineIter, warning=warning) : nothing

        if normΔs(robot) < ε && normf1 < ε
            foreach(s1tos0!,links)
            foreach(s1tos0!,constraints)
            return n
        else
            foreach(s1tos0!,links)
            foreach(s1tos0!,constraints)
            normf0=normf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf0,"."))
    end

    return -1
end

function lineSearch!(robot,normf0;iter=20, warning::Bool=false)
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
            return nothing
        end
    end

    if warning
        display(string("WARNING:  lineSearch! did not converge. n = ",iter,"."))
    end
end

lineStep!(node,diagonal,α) = (node.s1 = node.s0 - 1/(2^α)*diagonal.ŝ; nothing)
