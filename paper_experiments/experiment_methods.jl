function energy(link::Link,g,dt)
    x2 = link.x[2]
    v2 = 0.5*(getv1(link,dt)+getvnew(link))
    ω2 = 0.5*(getω1(link,dt)+getωnew(link))
    T = .5*(v2'*v2*link.m + ω2'*link.J*ω2)
    V = -link.m*g*(x2[3])
    return [T;V]
end

function energy(robot::Robot)
    en = [0.;0.]
    for link in robot.links
        en+=energy(link,robot.g,robot.dt)
    end
    return en
end

function drift(robot::Robot)
    d = 0
    for constraint in robot.constraints
        d += norm(g(constraint,robot))
    end
    return d
end

function simulate_energy!(robot::Robot;save::Bool=false,debug::Bool=false,disp::Bool=false)
    links = robot.links
    constraints = robot.constraints
    dt = robot.dt
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)

    totenergy = zeros(length(robot.steps),2)

    for i=robot.steps
        newton!(robot,warning=debug)
        save && saveToTraj!(robot,i)
        totenergy[i,:] = energy(robot)
        foreach(updatePos!,links,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return totenergy
end


function simulate_drift!(robot::Robot;save::Bool=false,debug::Bool=false,disp::Bool=false)
    links = robot.links
    constraints = robot.constraints
    dt = robot.dt
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)

    totdrift = zeros(length(robot.steps))

    for i=robot.steps
        newton!(robot,warning=debug)
        save && saveToTraj!(robot,i)
        totdrift[i] = drift(robot)
        foreach(updatePos!,links,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return totdrift
end

function simulate_reset!(robot::Robot,xin,qin;save::Bool=false,debug::Bool=false,disp::Bool=false,tol=1e-10)
    links = robot.links
    constraints = robot.constraints
    dt = robot.dt
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)

    for (i,link) in enumerate(links)
        link.x[1] = xin[i]
        link.x[2] = xin[i]
        link.q[1] = qin[i]
        link.q[2] = qin[i]
        link.s0 = @SVector zeros(6)
        link.s1 = @SVector zeros(6)
    end
    for constraint in constraints
        constraint.s0 *= 0.
        constraint.s1 *= 0.
    end

    for i=robot.steps
        newton!(robot,warning=debug,ε= tol)
        save && saveToTraj!(robot,i)
        foreach(updatePos!,links,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
end

function simulate_steptol!(robot::Robot;save::Bool=false,debug::Bool=false,disp::Bool=false)
    links = robot.links
    constraints = robot.constraints
    dt = robot.dt
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)

    stepandtol = zeros(100)

    for i=robot.steps
        stepandtol = newton_steptol!(robot,warning=debug)
        save && saveToTraj!(robot,i)
        foreach(updatePos!,links,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return stepandtol
end



function newton_steptol!(robot::Robot{T,Nl}; ε=1e-10, μ=1e-5, newtonIter=100, lineIter=10, warning::Bool=false) where {T,Nl}
    stepandtol = zeros(100)
    links = robot.links
    constraints = robot.constraints
    graph = robot.graph
    ldu = robot.ldu
    dt = robot.dt

    normf0 = normf(robot)
    for n=Base.OneTo(newtonIter)
        stepandtol[n]=normf0

        setentries!(robot)
        factor!(graph,ldu)
        solve!(graph,ldu) # x̂1 for each link and constraint

        foreach(update!,links,ldu)
        foreach(update!,constraints,ldu)
        # foreach(checkωnorm!,links,dt)

        normf1 = normf(robot)
        normf1>normf0 && lineSearch!(robot,normf0;iter=lineIter, warning=warning)



        if normΔs(robot) < ε && normf1 < ε
            stepandtol[n+1]=normf1
            foreach(s1tos0!,links)
            foreach(s1tos0!,constraints)
            return stepandtol
        else
            foreach(s1tos0!,links)
            foreach(s1tos0!,constraints)
            normf0=normf1
        end
    end

    if warning
        display(string("WARNING:  newton! did not converge. n = ",newtonIter,", tol = ",normf0,"."))
    end

    return stepandtol
end
