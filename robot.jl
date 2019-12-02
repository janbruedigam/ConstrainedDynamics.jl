using StaticArrays
using Base.Threads

include("util/quaternion.jl")
include("link.jl")
include("constraints/constraint.jl")
include("util/graph.jl")

mutable struct Robot{T,Nl,Nc,N}
    tend::T
    dt::T
    steps::UnitRange{Int64}

    nodes::Vector{Node{T}}
    nodesRange::Vector{UnitRange{Int64}}
    normf::T
    normDiff::T

    adjacency::Vector{SVector{N,Bool}}
    dfsgraph::Vector{SVector{N,Bool}}
    nodeList::SVector{N,Int64}
    parentList::SVector{N,Int64}
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T,Nl,Nc}) where {T,Nl,Nc}
    summary(io, R); println(io, " with ", Nl, " links and ", Nc, " constraints")
end


function Robot(links::Vector{<:Link{T}},constraints::Vector{<:Constraint{T}}; tend::T=10., dt::T=.01, g::T=9.81, root=1) where T
    Nl = length(links)
    Nc = length(constraints)
    N = Nl+Nc
    steps = Int(ceil(tend/dt))

    for (i,link) in enumerate(links)
        link.data.id = i
        link.dt = dt
        link.g = g
        link.trajectoryX = repeat([@SVector zeros(T,3)],steps)
        link.trajectoryQ = repeat([Quaternion{T}()],steps)
        link.trajectoryΦ = zeros(T,steps)
    end

    for (i,constraint) in enumerate(constraints)
        constraint.data.id = i+Nl
    end

    nodes = [links;constraints]
    nodesRange = [[1:Nl];[Nl+1:N]]
    normf = zero(T)
    normDiff = zero(T)

    adjacency = adjacencyMat(constraints,N)
    dfsgraph, list = dfs(adjacency,root=root)
    parentList = Vector{T}(undef,N)
    for i=1:N
        parentList[i] = parentNode(dfsgraph,i)
    end

    Robot{T,Nl,Nc,N}(tend,dt,1:steps,nodes,nodesRange,normf,normDiff,adjacency,dfsgraph,list,parentList)
end


function factor!(robot::Robot)
    list = robot.nodeList
    parentList = robot.parentList
    dfsgraph = robot.dfsgraph
    nodes = robot.nodes

    for n in list
        node = nodes[n]
        p = parentList[n]

        setD!(node)
        p!=0 ? setJ!(node,nodes[p]) : nothing
    end

    for n in list
        node = nodes[n]
        p = parentList[n]


        for (i,child) in enumerate(dfsgraph[n])
            child ? updateD!(node,nodes[i]) : nothing
        end
        invertD!(node)
        p!=0 ? updateJ!(node) : nothing
    end
end

function solve!(robot::Robot{T,Nl}) where {T,Nl}
    list = robot.nodeList
    parentList = robot.parentList
    dfsgraph = robot.dfsgraph
    adjacency = robot.adjacency
    nodes = robot.nodes
    nodesRange = robot.nodesRange

    for n in list
        node = nodes[n]

        setSol!(node)
        # (A) For extended equations
        # if n in robot.nodesRange[1]
        #     for (i,connected) in enumerate(adjacency[n])
        #         connected ? addGtλ!(node,nodes[i]) : nothing # this only changes sth for a link
        #     end
        # end
        # end (A)
        for (i,child) in enumerate(dfsgraph[n])
            child ? LSol!(node,nodes[i]) : nothing
        end
    end

    for n in reverse(list)
        node = nodes[n]
        p = parentList[n]

        DSol!(node)
        p!=0 ? USol!(node,nodes[p]) : nothing
    end
    # (A) For simplified equations
    for n = nodesRange[2]
        addλ0!(nodes[n])
    end
    # end (A)
end

function normf(robot::Robot{T}) where T
    robot.normf = 0
    nodes = robot.nodes

    foreach(setNormf!,nodes,robot)
    foreach(addNormf!,nodes,robot)
    return sqrt(robot.normf)
end

@inline addNormf!(node,robot::Robot) = (robot.normf += node.data.normf; nothing)

function normDiff(robot::Robot)
    robot.normDiff = 0
    nodes = robot.nodes

    foreach(setNormDiff!,nodes)
    foreach(addNormDiff!,nodes,robot)
    return sqrt(robot.normDiff)
end

@inline addNormDiff!(node,robot::Robot) = (robot.normDiff += node.data.normDiff; nothing)


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
    foreach(s0tos1!,robot.nodes)
    for i=robot.steps
        newton!(robot,warning=debug)
        for n=robot.nodesRange[1]
            save ? saveToTraj!(robot.nodes[n],i) : nothing
            updatePos!(robot.nodes[n])
        end

        disp && (i*robot.dt)%1<robot.dt*(1.0-.1) ? display(i*robot.dt) : nothing
    end
end

function trajSFunc(robot::Robot{T,Nl}) where {T,Nl}
    t = zeros(T,Nl,robot.steps[end])
    for i=1:Nl
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
