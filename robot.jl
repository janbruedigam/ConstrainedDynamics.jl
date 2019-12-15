mutable struct Robot{T,N,No}
    tend::T
    dt::T
    steps::Base.OneTo{Int64}

    origin::Origin{T}
    links::Vector{Link{T}}
    constraints::Vector{Constraint{T}}
    ldict::Dict{Int64,Int64}
    cdict::Dict{Int64,Int64}
    #???
    normf::T
    normΔs::T

    graph::Graph{N}

    # ldu::SparseLDU{}
    storage::Storage

    #TODO no constraints input
    function Robot(origin::Origin{T},links::Vector{Link{T}},constraints::Vector{<:Constraint{T}}; tend::T=10., dt::T=.01, g::T=-9.81, rootid=1, No=2) where T
        Nl = length(links)
        Nc = length(constraints)
        N = Nl+Nc
        steps = Int(ceil(tend/dt))

        ldict = Dict{Int64,Int64}()
        for (ind,link) in enumerate(links)
            push!(link.x, [link.x[1] for i=1:No-1]...)
            push!(link.q, [link.q[1] for i=1:No-1]...)
            push!(link.F, [link.F[1] for i=1:No-1]...)
            push!(link.τ, [link.τ[1] for i=1:No-1]...)

            ldict[link.id] = ind
        end

        cdict = Dict{Int64,Int64}()
        for (ind,constraint) in enumerate(constraints)
            cdict[constraint.id] = ind
        end

        resetGlobalID()

        normf = zero(T)
        normΔs = zero(T)

        graph = Graph(origin,links,constraints)
        # ldu = SparseLDU(graph)

        storage = Storage{T}(steps,Nl)

        new{T,N,No}(tend,dt,Base.OneTo(steps),origin,links,constraints,ldict,cdict,normf,normΔs,graph,storage)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.links), " links and ", length(R.constraints), " constraints")
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
    robot.nodes[robot.idlist[id]]
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

@inline function saveToTraj!(robot::Robot{T,N},t) where {T,N}
    No = robot.origin.No

    for i=robot.nodesrange[1]
        robot.storage.x[i][t]=robot.nodes[i].x[No]
        robot.storage.q[i][t]=robot.nodes[i].q[No]
    end
    return nothing
end

# @inline function saveToTraj!(robot::Robot,link::Link,i,t)
#     robot.storage.x[i][t]=link.x[2]
#     robot.storage.q[i][t]=link.q[2]
# end
#
# @inline function saveToTraj!(robot::Robot,constraint::Link,i,t)
#     robot.storage.λ[i][t]=constraint.data.s1
# end

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
        save && saveToTraj!(robot,i)
        for n=robot.nodesrange[1]
            updatePos!(nodes[n])
        end

        disp && (i*robot.dt)%1<robot.dt*(1.0-.1) && display(i*robot.dt)
    end
end

function plotTraj(robot,trajS,id)
    p = plot(collect(0:robot.dt:robot.tend-robot.dt),trajS[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:robot.dt:robot.tend-robot.dt),trajS[ind,:])
    end
    return p
end
