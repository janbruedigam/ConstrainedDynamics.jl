mutable struct Robot{T,N,NpF}
    tend::T
    dt::T
    steps::UnitRange{Int64}

    origin::Link{T,0,0}
    nodes::Vector{Node{T}}
    fillins::Vector{FillIn{T}}
    nodesrange::Vector{UnitRange{Int64}}
    normf::T
    normΔs::T

    idlist::SVector{NpF,Int64}

    graph::Graph{N}

    function Robot(origin::Link{T,0,0},links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}}; tend::T=10., dt::T=.01, g::T=-9.81, rootid=1) where T
        Nl = length(links)
        Nc = length(constraints)
        N = Nl+Nc+1
        steps = Int(ceil(tend/dt))

        origin.dt = dt
        origin.g = g
        origin.data.id = 1
        for (i,link) in enumerate(links)
            link.data.id = i+1
            link.dt = dt
            link.g = g
            link.trajectoryX = repeat([@SVector zeros(T,3)],steps)
            link.trajectoryQ = repeat([Quaternion{T}()],steps)
            link.trajectoryΦ = zeros(T,steps)
        end

        for (i,constraint) in enumerate(constraints)
            constraint.data.id = i+1+Nl
        end

        #TODO do constraint properly
        resetGlobalID()

        nodes = [links;constraints]

        nodesrange = [[1:Nl];[Nl+1:Nl+Nc]]
        normf = zero(T)
        normΔs = zero(T)

        graph = Graph(origin,links,constraints)

        idlist = zeros(Int64,N)
        idlist[origin.data.id] = 0
        for i=1:N-1
            idlist[nodes[i].data.id] = i
        end

        fillins = createfillins(graph,idlist,origin,nodes)
        F = length(fillins)

        idlist = [idlist;zeros(Int64,F)]

        for i=1:F
            idlist[fillins[i].data.id] = i
        end

        new{T,N,N+F}(tend,dt,1:steps,origin,nodes,fillins,nodesrange,normf,normΔs,idlist,graph)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.nodesrange[1]), " links and ", length(R.nodesrange[2]), " constraints")
end


function factor!(robot::Robot)
    graph = robot.graph
    list = graph.dfslist
    pattern = graph.pattern

    for id in list
        if !isroot(graph,id)
            node = getnode(robot,id)

            for cid in list # in correct order
                cid==id && break
                if pattern[id][cid]!=0 # is actually a (loop) child
                    cfillin = getfillin(robot,cid,id)
                    childnode = getnode(robot,cid)
                    setJ!(childnode,node,cfillin)
                    for gcid in list # in correct order
                        gcid==cid && break
                        if pattern[id][gcid]!=0 && pattern[cid][gcid]!=0# is actually a (loop child)
                            gcfillin = getfillin(robot,gcid,id)
                            cgcfillin = getfillin(robot,gcid,cid)
                            grandchildnode = getnode(robot,gcid)
                            updateJ1!(grandchildnode,gcfillin,cgcfillin,cfillin)
                        end
                    end
                    updateJ2!(childnode,cfillin)
                end
            end

            setD!(node)
            for (cid,ischild) in enumpattern(graph,id)
                ischild!=0 && updateD!(node,getnode(robot,cid),getfillin(robot,cid,id))
            end
            invertD!(node)
        end
    end
end

function solve!(robot::Robot)
    graph = robot.graph
    list = graph.dfslist
    nodes = robot.nodes
    pattern = graph.pattern

    for id in list
        if !isroot(graph,id)
            node = getnode(robot,id)

            setSol!(node)

            for cid in list # in correct order (shouldnt matter here)
                cid==id && break
                if pattern[id][cid]!=0 # is actually a (loop) child
                    LSol!(node,getnode(robot,cid),getfillin(robot,cid,id))
                end
            end
        end
    end

    for id in reverse(list)
        if !isroot(graph,id)
            node = getnode(robot,id)

            DSol!(node)

            for pid in reverse(list)
                pid==id && break
                if pattern[pid][id]!=0 # is a (loop) parent
                    USol!(node,getnode(robot,pid),getfillin(robot,id,pid))
                end
            end
        end
    end

    # (A) For simplified equations
    for indn = robot.nodesrange[2]
        addλ0!(nodes[indn])
    end
    # end (A)
end

@inline function getnode(robot::Robot,id::Int64)
    isroot(robot.graph,id) ? robot.origin : robot.nodes[robot.idlist[id]]
end

@inline function getfillin(robot::Robot,cid::Int64,pid::Int64)
    fid = robot.graph.pattern[pid][cid]
    robot.fillins[robot.idlist[fid]]
end


function normf(robot::Robot{T}) where T
    robot.normf = 0
    nodes = robot.nodes

    foreach(setNormf!,nodes,robot)
    foreach(addNormf!,nodes,robot)
    return sqrt(robot.normf)
end

@inline addNormf!(node,robot::Robot) = (robot.normf += node.data.normf; nothing)

function normΔs(robot::Robot)
    robot.normΔs = 0
    nodes = robot.nodes

    foreach(setNormΔs!,nodes)
    foreach(addNormΔs!,nodes,robot)
    return sqrt(robot.normΔs)
end

@inline addNormΔs!(node,robot::Robot) = (robot.normΔs += node.data.normΔs; nothing)


@inline function saveToTraj!(link::Link,i)
    No = link.No
    link.trajectoryX[i] = link.x[No]
    link.trajectoryQ[i] = link.q[No]
    link.trajectoryΦ[i] = angleAxis(link.q[No])[1]*sign(angleAxis(link.q[No])[2][1])
    return nothing
end

@inline function updatePos!(link::Link)
    link.x[1] = link.x[2]
    link.x[2] += getvnew(link)*link.dt
    link.q[1] = link.q[2]
    link.q[2] = link.dt/2*(Lmat(link.q[2])*ωbar(link))
    return nothing
end


function sim!(robot::Robot;save::Bool=false,debug::Bool=false,disp::Bool=false)
    nodes = robot.nodes
    foreach(s0tos1!,nodes)
    for i=robot.steps
        newton!(robot,warning=debug)
        for n=robot.nodesrange[1]
            save ? saveToTraj!(nodes[n],i) : nothing
            updatePos!(nodes[n])
        end

        disp && (i*robot.dt)%1<robot.dt*(1.0-.1) ? display(i*robot.dt) : nothing
    end
end

function trajSFunc(robot::Robot{T}) where {T}
    t = zeros(T,length(robot.nodesrange[1]),length(robot.steps))
    for i=robot.nodesrange[1]
            t[i,:] = robot.nodes[i].trajectoryΦ
    end
    return t
end

function plotTraj(robot,trajS,id)
    p = plot(collect(0:robot.dt:robot.tend-robot.dt),trajS[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:robot.dt:robot.tend-robot.dt),trajS[ind,:])
    end
    return p
end
