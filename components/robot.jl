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

        storage = Storage{T}(steps,Nl,Nc)

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
            setJ!(robot,getentry(ldu,(id,cid)),id,getconstraint(robot,cid))
        end

        diagonal = getentry(ldu,id)
        link = getlink(robot,id)
        setDandŝ!(diagonal,link,robot)
    end

    for node in robot.constraints
        id = node.id

        for cid in directchildren(graph,id)
            setJ!(robot,getentry(ldu,(id,cid)),node,cid)
        end

        for cid in loopchildren(graph,id)
            setJ!(getentry(ldu,(id,cid)))
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,node,robot)
    end
end

@inline function correctλ!(robot::Robot)
    for constraint in robot.constraints
        addλ0!(getentry(robot.ldu,constraint.id),constraint)
    end
end

@inline getlink(robot::Robot,id::Int64) = robot.links[robot.ldict[id]]
@inline getlink(robot::Robot,id::Nothing) = robot.origin
@inline getconstraint(robot::Robot,id::Int64) = robot.constraints[robot.cdict[id]]
# @inline function getnode(robot::Robot,id::Int64) # should only be used in setup
#      if haskey(robot.ldict,id)
#          return getlink(robot,id)
#      elseif haskey(robot.cdict,id)
#          return getconstraint(robot,id)
#      elseif id == robot.originid
#          return robot.origin
#      else
#          error("not found.")
#      end
#  end


@inline function normf(robot::Robot{T}) where T
    robot.normf = 0

    for link in robot.links
        robot.normf+=normf(link,robot)
    end
    foreach(addNormf!,robot.constraints,robot)

    return sqrt(robot.normf)
end

@inline addNormf!(node,robot::Robot) = (robot.normf += normf(node,robot); nothing)

@inline function normΔs(robot::Robot)
    robot.normΔs = 0

    robot.normΔs+=mapreduce(normΔs,+,robot.links)
    foreach(addNormΔs!,robot.constraints,robot)

    return sqrt(robot.normΔs)
end

@inline addNormΔs!(node,robot::Robot) = (robot.normΔs += normΔs(node); return)

function saveToTraj!(robot::Robot,t)
    No = robot.No
    for (ind,link) in enumerate(robot.links)
        robot.storage.x[ind][t]=link.x[No]
        robot.storage.q[ind][t]=link.q[No]
    end
    for (ind,constraint) in enumerate(robot.constraints)
        robot.storage.λ[ind][t]=constraint.s1
    end
    return nothing
end

@inline function updatePos!(link::Link,dt)
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

function plotθ(robot::Robot{T},id) where T
    n = length(robot.links)
    θ = zeros(T,n,length(robot.steps))
    for i=1:n
        qs = robot.storage.q[i]
        for (t,q) in enumerate(qs)
            θ[i,t] = angleAxis(q)[1]*sign(angleAxis(q)[2][1])
        end
    end

    p = plot(collect(0:robot.dt:robot.tend-robot.dt),θ[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:robot.dt:robot.tend-robot.dt),θ[ind,:])
    end
    return p
end

function plotλ(robot::Robot{T},id) where T
    n = sum(length.(robot.constraints))
    λ = zeros(T,n,length(robot.steps))
    startpos = 1
    endpos = 0
    for i=1:length(robot.constraints)
        endpos = startpos + length(robot.constraints[i]) -1

        λs = robot.storage.λ[i]
        for (t,val) in enumerate(λs)
            λ[startpos:endpos,t] = val
        end

        startpos = endpos + 1
    end

    p = plot(collect(0:robot.dt:robot.tend-robot.dt),λ[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:robot.dt:robot.tend-robot.dt),λ[ind,:])
    end
    return p
end
