mutable struct Robot{T,N,No}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T

    origin::Link{T,0}
    nodes::Vector{Node{T}}
    nodesrange::Vector{UnitRange{Int64}}
    dict::Dict{Int64,Int64}

    #???
    normf::T
    normΔs::T

    graph::Graph{N}

    ldu::SparseLDU{T}
    storage::Storage{T}

    #TODO no constraints input
    function Robot(origin::Link{T,0},links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}}; tend::T=10., dt::T=.01, g::T=-9.81, rootid=1, No=2) where T
        Nl = length(links)
        Nc = length(constraints)
        N = Nl+Nc
        steps = Int(ceil(tend/dt))

        ldict = Dict{Int64,Int64}()

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

            ldict[link.id] = ind
        end

        cdict = Dict{Int64,Int64}()
        for (ind,constraint) in enumerate(constraints)
            cdict[constraint.id] = ind
        end

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
        ldu = SparseLDU(graph,links,constraints,ldict,cdict)

        storage = Storage{T}(steps,Nl)

        new{T,N,No}(tend,Base.OneTo(steps),dt,g,origin,nodes,nodesrange,dict,normf,normΔs,graph,ldu,storage)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.nodesrange[1]), " links and ", length(R.nodesrange[2]), " constraints")
end

function setentries!(robot::Robot)
    ldu = robot.ldu
    graph = robot.graph

    for id in graph.dfslist
        node = getnode(robot,id)
        diagonal = getentry(ldu,id)

        for cid in successors(graph,id)
            cid == -1 && break
            offdiagonal = getentry(ldu,(id,cid))
            cnode =  getnode(robot,cid)
            setJ!(offdiagonal,node,cnode)
        end

        setD!(diagonal,node)
        setSol!(diagonal,node)
    end
end

@inline function correctλ!(robot::Robot)
    nodes = robot.nodes
    diagonals = robot.ldu.diagonals
    for indn = robot.nodesrange[2]
        addλ0!(diagonals[indn],nodes[indn])
    end
end

@inline getnode(robot::Robot,id::Int64) = robot.nodes[robot.dict[id]]

function normf(robot::Robot{T}) where T
    robot.normf = 0
    nodes = robot.nodes

    foreach(setNormf!,nodes,robot)
    foreach(addNormf!,nodes,robot)
    return sqrt(robot.normf)
end

@inline addNormf!(node,robot::Robot) = (robot.normf += node.data.normf; nothing)

@inline function normΔs(robot::Robot)
    robot.normΔs = 0
    nodes = robot.nodes

    foreach(setNormΔs!,nodes)
    foreach(addNormΔs!,nodes,robot)
    return sqrt(robot.normΔs)
end

@inline addNormΔs!(node,robot::Robot) = (robot.normΔs += node.data.normΔs; nothing)

@inline function saveToTraj!(robot::Robot{T,N,No},t) where {T,N,No}
    for i=robot.nodesrange[1]
        robot.storage.x[i][t]=robot.nodes[i].x[No]
        robot.storage.q[i][t]=robot.nodes[i].q[No]
    end
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
