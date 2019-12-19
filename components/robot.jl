mutable struct Robot{T,N}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T
    No::Int64

    origin::Origin{T}
    links::Vector{Link{T}}
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
        ldu = SparseLDU(graph,links,constraints,ldict,cdict)

        storage = Storage{T}(steps,Nl)

        new{T,N}(tend,Base.OneTo(steps),dt,g,No,origin,links,constraints,ldict,cdict,normf,normΔs,graph,ldu,storage)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.links), " links and ", length(R.constraints), " constraints")
end

function setentries!(robot::Robot)

    graph = robot.graph
    ldu = robot.ldu

    for (id,_) in robot.ldict
        for cid in directchildren(graph,id)
            # cid == -1 && break
            setJ!(robot,getentry(ldu,(id,cid)),id,getconstraint(robot,cid))
        end

        diagonal = getentry(ldu,id)
        link = getlink(robot,id)
        setDandŝ!(diagonal,link,robot)
    end

    for node in robot.constraints
        id = node.id

        for cid in directchildren(graph,id)
            # cid == -1 && break
            setJ!(robot,getentry(ldu,(id,cid)),node,cid)
        end

        for cid in loopchildren(graph,id)
            # cid == -1 && break
            setJ!(getentry(ldu,(id,cid)))
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,node,robot)
    end
end

function correctλ!(robot::Robot)
    for constraint in robot.constraints
        addλ0!(getentry(robot.ldu,constraint.id),constraint)
    end
end

getlink(robot::Robot,id::Int64) = robot.links[robot.ldict[id]]
getlink(robot::Robot,id::Nothing) = robot.origin
getconstraint(robot::Robot,id::Int64) = robot.constraints[robot.cdict[id]]


function normf(robot::Robot{T}) where T
    robot.normf = 0

    for link in robot.links
        robot.normf+=normf(link,robot)
    end
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

function saveToTraj!(robot::Robot,t)
    No = robot.No
    for (ind,link) in enumerate(robot.links)
        robot.storage.x[ind][t]=link.x[No]
        robot.storage.q[ind][t]=link.q[No]
    end
    return nothing
end

function updatePos!(link::Link,dt)
    link.x[1] = link.x[2]
    link.x[2] += getvnew(link)*dt
    link.q[1] = link.q[2]
    link.q[2] = dt/2*(Lmat(link.q[2])*ωbar(link,dt))
    return nothing
end


function simulate!(robot::Robot;save::Bool=false,debug::Bool=false,disp::Bool=false)
    links = robot.links
    constraints = robot.constraints
    dt = robot.dt
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)
    for i=robot.steps
        newton!(robot,warning=debug)
        save && saveToTraj!(robot,i)
        for link in links
            updatePos!(link,dt)
        end

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
end

function plotTraj(robot::Robot{T},id) where T
    n = length(robot.links)
    angles = zeros(T,n,length(robot.steps))
    for i=1:n
        qs = robot.storage.q[i]
        for (t,q) in enumerate(qs)
            angles[i,t] = angleAxis(q)[1]*sign(angleAxis(q)[2][1])
        end
    end

    p = plot(collect(0:robot.dt:robot.tend-robot.dt),angles[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:robot.dt:robot.tend-robot.dt),angles[ind,:])
    end
    return p
end
