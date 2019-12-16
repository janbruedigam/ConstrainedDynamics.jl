function newton!(robot::Robot{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=20, warning::Bool=false) where {T,Nl}
    n = 1
    nodes = robot.nodes
    diagonals = robot.diagonals
    normf0 = normf(robot)
    for outer n=1:newtonIter
        setentries!(robot)
        factor!(robot)
        solve!(robot) # x̂1 for each link and constraint
        # foreach(update!,nodes) # x1 = x0 - x̂1 for each link and constraint
        for (i,node) in enumerate(nodes)
            update!(node,diagonals[i])
        end

        normf1 = normf(robot)
        normf1>normf0 ? lineSearch!(robot,normf0;iter=lineIter, warning=warning) : nothing

        if normΔs(robot) < ε && normf1 < ε
            foreach(s1tos0!,nodes)
            return n
        else
            foreach(s1tos0!,nodes)
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
    nodes = robot.nodes
    diagonals = robot.diagonals
    for (i,node) in enumerate(nodes)
        lineStep!(node,diagonals[i],α)# x1 = x0 + 1/(2^α)*d
    end

    for n=1:iter
        α += 1
        if normf(robot) >= normf0
            for (i,node) in enumerate(nodes)
                lineStep!(node,diagonals[i],α)# x1 = x0 + 1/(2^α)*d
            end
        else
            return nothing
        end
    end

    if warning
        display(string("WARNING:  lineSearch! did not converge. n = ",iter,"."))
    end
end

@inline lineStep!(node,diagonal,α) = (d = node.data; d.s1 = d.s0 - 1/(2^α)*diagonal.ŝ; nothing)
