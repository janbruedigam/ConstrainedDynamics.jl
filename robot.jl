mutable struct Robot{T,N,No}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T

    origin::Link{T,0,0}
    # links::Vector{Link{T}}
    # constraints::Vector{Constraint{T}}
    nodes::Vector{Node{T}}
    fillins::Vector{OffDiagonalEntry{T}}
    nodesrange::Vector{UnitRange{Int64}}
    # ldict::Dict{Int64,Int64}
    # cdict::Dict{Int64,Int64}
    dict::Dict{Int64,Int64}
    fdict::Dict{Tuple{Int64,Int64},Int64}
    #???
    normf::T
    normΔs::T

    graph::Graph{N}

    # ldu::SparseLDU{T}
    storage::Storage{T}

    #TODO no constraints input
    function Robot(origin::Link{T,0,0},links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}}; tend::T=10., dt::T=.01, g::T=-9.81, rootid=1, No=2) where T
        Nl = length(links)
        Nc = length(constraints)
        N = Nl+Nc
        steps = Int(ceil(tend/dt))

        # ldict = Dict{Int64,Int64}()

        origin.g = g
        origin.dt = dt
        origin.No = No
        push!(origin.x, [origin.x[1] for i=1:No-1]...)
        push!(origin.q, [origin.q[1] for i=1:No-1]...)
        push!(origin.F, [origin.F[1] for i=1:No-1]...)
        push!(origin.τ, [origin.τ[1] for i=1:No-1]...)
        for (ind,link) in enumerate(links)
            link.g = g
            link.dt = dt
            link.No = No
            push!(link.x, [link.x[1] for i=1:No-1]...)
            push!(link.q, [link.q[1] for i=1:No-1]...)
            push!(link.F, [link.F[1] for i=1:No-1]...)
            push!(link.τ, [link.τ[1] for i=1:No-1]...)

            # ldict[link.id] = ind
        end

        # cdict = Dict{Int64,Int64}()
        # for (ind,constraint) in enumerate(constraints)
        #     cdict[constraint.id] = ind
        # end

        nodes = [links;constraints]
        nodesrange = [[1:Nl];[Nl+1:Nl+Nc]]
        dict = Dict{Int64,Int64}()
        for (ind,node) in enumerate(nodes)
            dict[node.id] = ind
        end

        resetGlobalID()

        normf = zero(T)
        normΔs = zero(T)

        graph = Graph(origin,links,constraints)
        # ldu = SparseLDU(graph,links,constraints,ldict,cdict)

        fillins,fdict = createfillins(graph,dict,origin,nodes)

        storage = Storage{T}(steps,Nl)

        new{T,N,No}(tend,Base.OneTo(steps),dt,g,origin,nodes,fillins,nodesrange,dict,fdict,normf,normΔs,graph,storage)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.nodesrange[1]), " links and ", length(R.nodesrange[2]), " constraints")
end

# function setentries!(robot::Robot{T,N,No}) where {T,N,No}
#     ldu = robot.ldu
#     ldict = robot.ldict
#     cdict = robot.cdict
#
#     links = robot.links
#     constraints = robot.constraints
#     dt = robot.dt
#     g = robot.g
#
#     graph = robot.graph
#     dfsgraph = graph.dfsgraph
#     fillins = graph.fillins
#     dict = graph.dict
#     rdict = graph.rdict
#
#     diagonals = ldu.diagonals
#     offdiagonals = ldu.offdiagonals
#
#     for link in links
#         id = link.id
#         entry = getentry(ldu,id)
#         setD!(entry,link,dt)
#         # setSol!(entry,link,dt,g,No)
#         # for (ind,ischild) in enumerate(dfsgraph[dict[id]])
#         #     if ischild
#         #         childid = rdict[ind]
#         #         childind = cdict[childid]
#         #         setJ!(getentry(ldu,(id,childid)),link,constraints[childind],dt,No)
#         #     end
#         # end
#         # for (ind,ischild) in enumerate(fillins[dict[id]])
#         #     if ischild
#         #         childid = rdict[ind]
#         #         setJ!(getentry(ldu,(id,childid)))
#         #     end
#         # end
#     end
#
#     # for constraint in constraints
#     #     id = constraint.id
#     #     entry = getentry(ldu,id)
#     #     setD!(entry,constraint)
#     #     setSol!(entry,constraint,dt)
#     #     for (ind,ischild) in enumerate(dfsgraph[dict[id]])
#     #         if ischild
#     #             childid = rdict[ind]
#     #             childind = ldict[childid]
#     #             setJ!(getentry(ldu,(id,childid)),constraint,links[childind],dt,No)
#     #         end
#     #     end
#     #     for (ind,ischild) in enumerate(fillins[dict[id]])
#     #         if ischild
#     #             childid = rdict[ind]
#     #             setJ!(getentry(ldu,(id,childid)))
#     #         end
#     #     end
#     # end
# end

function createfillins(graph::Graph,dict,origin::Link{T},nodes::Vector) where T
    rdict = graph.rdict
    fdict = Dict{Tuple{Int64,Int64},Int64}()
    fillins = Vector{OffDiagonalEntry}(undef,0)

    indcounter = 0

    for (parentind,row) in enumerate(graph.pattern)
        for (childind,ischild) in enumerate(row)
            if ischild
                parentnode = nodes[dict[rdict[parentind]]]
                childnode = nodes[dict[rdict[childind]]]

                fdict[(rdict[parentind],rdict[childind])] = indcounter+=1

                push!(fillins,OffDiagonalEntry{T,length(childnode),length(parentnode)}())
            end
        end
    end
    return fillins,fdict
end

function factor!(robot::Robot)
    graph = robot.graph
    list = graph.dfslist
    pattern = graph.pattern
    nodes = robot.nodes
    dict = robot.dict
    gdict = graph.dict
    grdict = graph.rdict
    fdict = robot.fdict
    fillins = robot.fillins

    for id in list
        node = nodes[dict[id]]#getnode(robot,id)

        for cid in list # in correct order
            cid==id && break
            if pattern[gdict[id]][gdict[cid]] # is actually a (loop) child
                cfillin = fillins[fdict[(id,cid)]]#getfillin(robot,cid,id)
                childnode = nodes[dict[cid]]#getnode(robot,cid)
                setJ!(childnode,node,cfillin)
                for gcid in list # in correct order
                    gcid==cid && break
                    if pattern[gdict[id]][gdict[gcid]] && pattern[gdict[cid]][gdict[gcid]] # is actually a (loop child)
                        gcfillin = fillins[fdict[(id,gcid)]]#getfillin(robot,gcid,id)
                        cgcfillin = fillins[fdict[(cid,gcid)]]#getfillin(robot,gcid,cid)
                        grandchildnode = nodes[dict[gcid]]#getnode(robot,gcid)
                        updateJ1!(grandchildnode,gcfillin,cgcfillin,cfillin)
                    end
                end
                updateJ2!(childnode,cfillin)
            end
        end

        setD!(node)
        for (cind,ischild) in enumerate(pattern[gdict[id]])
            # ischild!=0 && updateD!(node,getnode(robot,cid),getfillin(robot,cid,id))
            if ischild
                cid = grdict[cind]
                updateD!(node,nodes[dict[cid]],fillins[fdict[(id,cid)]])
            end
        end
        invertD!(node)
    end
end

function solve!(robot::Robot)
    graph = robot.graph
    list = graph.dfslist
    nodes = robot.nodes
    pattern = graph.pattern
    nodes = robot.nodes
    dict = robot.dict
    fdict = robot.fdict
    gdict = graph.dict
    grdict = graph.rdict
    fillins = robot.fillins

    for id in list
        node = nodes[dict[id]]#getnode(robot,id)

        setSol!(node)

        for cid in list # in correct order (shouldnt matter here)
            cid==id && break
            if pattern[gdict[id]][gdict[cid]] # is actually a (loop) child
                LSol!(node,nodes[dict[cid]],fillins[fdict[(id,cid)]])
            end
        end
    end

    for id in reverse(list)
        node = nodes[dict[id]]#getnode(robot,id)

        DSol!(node)

        for pid in reverse(list)
            pid==id && break
            if pattern[gdict[pid]][gdict[id]] # is a (loop) parent
                USol!(node,nodes[dict[pid]],fillins[fdict[(pid,id)]])
            end
        end
    end

    # (A) For simplified equations
    for indn = robot.nodesrange[2]
        addλ0!(nodes[indn])
    end
    # end (A)
end

# @inline function getnode(robot::Robot,id::Int64)
#     robot.nodes[robot.idlist[id]]
# end
#
# @inline function getfillin(robot::Robot,cid::Int64,pid::Int64)
#     fid = robot.graph.pattern[pid][cid]
#     robot.fillins[robot.idlist[fid]]
# end


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
