mutable struct Robot{T,N}
    tend::T
    dt::T
    steps::UnitRange{Int64}

    origin::Link{T,0,0,0,0}
    nodes::Vector{Node{T}}
    nodesrange::Vector{UnitRange{Int64}}
    normf::T
    normΔs::T

    idlist::SVector{N,Int64} # matches ids to "nodes" index, 0 is for root

    graph::Graph{N}
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(robot.nodesrange[1]), " links and ", length(robot.nodesrange[2]), " constraints")
end


function Robot(origin::Link{T,0,0,0,0},links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}}; tend::T=10., dt::T=.01, g::T=-9.81, rootid=1) where T
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

    nodes = [links;constraints]

    idlist = zeros(Int64,N)
    idlist[1] = 0
    for i=1:N-1
        idlist[i+1] = i
    end

    nodesrange = [[1:Nl];[Nl+1:Nl+Nc]]
    normf = zero(T)
    normΔs = zero(T)

    graph = Graph(constraints)

    Robot{T,N}(tend,dt,1:steps,origin,nodes,nodesrange,normf,normΔs,idlist,graph)
end


function factor!(robot::Robot)
    graph = robot.graph
    list = graph.dfslist

    for id in list
        if !isroot(graph,id)
            node = getnode(robot,id)
            pid = parent(graph,id)

            setD!(node)
            !isroot(graph,pid) && setJ!(node,getnode(robot,pid))
        end
    end

    for id in list
        if !isroot(graph,id)
            node = getnode(robot,id)

            for (cid,ischild) in enumchildren(graph,id)
                ischild && updateD!(node,getnode(robot,cid))
            end
            invertD!(node)
            updateJ!(node)
        end
    end
end

function solve!(robot::Robot)
    graph = robot.graph
    list = graph.dfslist
    nodes = robot.nodes

    for id in list
        if !isroot(graph,id)
            node = getnode(robot,id)

            setSol!(node)
            for (cid,ischild) in enumchildren(graph,id)
                ischild && LSol!(node,getnode(robot,cid))
            end
        end
    end

    for id in reverse(list)
        if !isroot(graph,id)
            node = getnode(robot,id)
            pid = parent(graph,id)

            DSol!(node)
            !isroot(graph,pid) && USol!(node,getnode(robot,pid))
        end
    end

    # (A) For simplified equations
    for indn = robot.nodesrange[2]
        addλ0!(nodes[indn])
    end
    # end (A)
end

@inline function getnode(robot::Robot,id::Int64)
    graph = robot.graph
    isroot(graph,id) ? robot.origin : robot.nodes[robot.idlist[id]]
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
