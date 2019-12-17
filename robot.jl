mutable struct Robot{T,N,No}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T

    origin::Link{T,0}
    links::Vector{Link{T,6}}
    constraints::Vector{Constraint{T}}
    ldict::Dict{Int64,Int64}
    cdict::Dict{Int64,Int64}

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

        resetGlobalID()

        normf = zero(T)
        normΔs = zero(T)

        graph = Graph(origin,links,constraints)
        ldu = SparseLDU(graph,links,constraints,ldict,cdict)

        storage = Storage{T}(steps,Nl)

        new{T,N,No}(tend,Base.OneTo(steps),dt,g,origin,links,constraints,ldict,cdict,normf,normΔs,graph,ldu,storage)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.links[1]), " links and ", length(R.constraints[2]), " constraints")
end

function setentries!(robot::Robot)
    setentries!(robot,robot.links)
    setentries!(robot,robot.constraints)
end

function correctλ!(robot::Robot)
    for constraint in robot.constraints
        addλ0!(getentry(robot.ldu,constraint.id),constraint)
    end
end

getlink(robot::Robot,id::Int64) = robot.links[robot.ldict[id]]
getconstraint(robot::Robot,id::Int64) = robot.constraints[robot.cdict[id]]
getnode(robot::Robot,id::Int64) = haskey(robot.ldict,id) ? getlink(robot,id) : getconstraint(robot,id)

function normf(robot::Robot{T}) where T
    robot.normf = 0

    foreach(addNormf!,robot.links,robot)
    foreach(addNormf!,robot.constraints,robot)
    
    return sqrt(robot.normf)
end

addNormf!(node,robot::Robot) = (robot.normf += normf(node,robot); nothing)

function normΔs(robot::Robot)
    robot.normΔs = 0

    robot.normΔs+=mapreduce(normΔs,+,robot.links)
    foreach(addNormΔs!,robot.constraints,robot)

    return sqrt(robot.normΔs)
end

addNormΔs!(node,robot::Robot) = (robot.normΔs += normΔs(node); return)

function saveToTraj!(robot::Robot{T,N,No},t) where {T,N,No}
    for (ind,link) in enumerate(robot.links)
        robot.storage.x[ind][t]=link.x[No]
        robot.storage.q[ind][t]=link.q[No]
    end
    return nothing
end

function updatePos!(link::Link)
    link.x[1] = link.x[2]
    link.x[2] += getvnew(link)*link.dt
    link.q[1] = link.q[2]
    link.q[2] = link.dt/2*(Lmat(link.q[2])*ωbar(link))
    return nothing
end


function sim!(robot::Robot;save::Bool=false,debug::Bool=false,disp::Bool=false)
    links = robot.links
    constraints = robot.constraints
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)
    for i=robot.steps
        newton!(robot,warning=debug)
        save && saveToTraj!(robot,i)
        for link in links
            updatePos!(link)
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
