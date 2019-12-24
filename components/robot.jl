mutable struct Robot{T,N}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T
    No::Int64

    origin::Origin{T}
    links::UnitDict{Base.OneTo{Int64},Link{T}}
    constraints::UnitDict{UnitRange{Int64},<:Constraint{T}}

    #TODO remove once Constraint is homogenous
    normf::T
    normΔs::T

    graph::Graph{N}

    ldu::SparseLDU{T}
    storage::Storage{T}

    #TODO no constraints input
    function Robot(origin::Origin{T},links::Vector{Link{T}},constraints::Vector{<:Constraint{T}};
        tend::T=10., dt::T=.01, g::T=-9.81, No=2) where T


        resetGlobalID()

        Nl = length(links)
        Nc = length(constraints)
        N = Nl+Nc
        steps = Int(ceil(tend/dt))

        currentid = 1

        ldict = Dict{Int64,Int64}()
        for (ind,link) in enumerate(links)
            push!(link.x, [link.x[1] for i=1:No-1]...)
            push!(link.q, [link.q[1] for i=1:No-1]...)
            push!(link.F, [link.F[1] for i=1:No-1]...)
            push!(link.τ, [link.τ[1] for i=1:No-1]...)

            for c in constraints
                c.pid == link.id && (c.pid = currentid)
                for (ind,linkid) in enumerate(c.linkids)
                    if linkid == link.id
                        c.linkids = setindex(c.linkids,currentid,ind)
                        c.constraints[ind].cid = currentid
                    end
                end
            end

            link.id = currentid
            currentid+=1

            ldict[link.id] = ind
        end

        cdict = Dict{Int64,Int64}()
        for (ind,constraint) in enumerate(constraints)
            constraint.id = currentid
            currentid+=1

            cdict[constraint.id] = ind
        end

        normf = zero(T)
        normΔs = zero(T)

        graph = Graph(origin,links,constraints)
        ldu = SparseLDU(graph,links,constraints,ldict,cdict)

        storage = Storage{T}(steps,Nl,Nc)

        links = UnitDict(links)
        constraints = UnitDict((links[Nl].id+1):currentid-1,constraints)

        new{T,N}(tend,Base.OneTo(steps),dt,g,No,origin,links,constraints,normf,normΔs,graph,ldu,storage)
    end

    function Robot(origin::Origin{T},links::Vector{Link{T}};
        tend::T=10., dt::T=.01, g::T=-9.81, No=2) where T

        constraints = Vector{Constraint{T}}(undef,0)
        for link in links
            push!(constraints,Constraint(OriginConnection(origin,link)))
        end
        Robot(origin,links,constraints,tend=tend, dt=dt, g=g, No=No)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Robot{T}) where {T}
    summary(io, R); println(io, " with ", length(R.links), " links and ", length(R.constraints), " constraints")
end

function setentries!(robot::Robot)
    graph = robot.graph
    ldu = robot.ldu

    for (id,link) in pairs(robot.links)
        for cid in directchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)),id,getconstraint(robot,cid),robot)
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,link,robot)
    end

    for node in robot.constraints
        id = node.id

        for cid in directchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)),node,cid,robot)
        end

        for cid in loopchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)))
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,node,robot)
    end
end

@inline getlink(robot::Robot,id::Int64) = robot.links[id]
@inline getlink(robot::Robot,id::Nothing) = robot.origin
@inline getconstraint(robot::Robot,id::Int64) = robot.constraints[id]

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

@inline function normf(link::Link{T},robot::Robot) where T
    f = dynamics(link,robot)
    return dot(f,f)
end

@inline function normf(c::Constraint,robot::Robot)
    f = g(c,robot)
    return dot(f,f)
end

@inline function GtλTof!(link::Link,c::Constraint,robot)
    link.f -= ∂g∂pos(c,link.id,robot)'*c.s1
    return
end

@inline function normf(robot::Robot)
    robot.normf = 0

    for link in robot.links
        robot.normf+=normf(link,robot)
    end
    foreach(addNormf!,robot.constraints,robot)

    return sqrt(robot.normf)
end

@inline function normΔs(robot::Robot)
    robot.normΔs = 0

    robot.normΔs += mapreduce(normΔs,+,robot.links)
    foreach(addNormΔs!,robot.constraints,robot)

    return sqrt(robot.normΔs)
end

@inline function addNormf!(c::Constraint,robot::Robot)
    robot.normf += normf(c,robot)
    return
end

@inline function addNormΔs!(component::Component,robot::Robot)
    robot.normΔs += normΔs(component)
    return
end

function saveToTraj!(robot::Robot,t)
    No = robot.No
    for (ind,link) in enumerate(robot.links)
        robot.storage.x[ind][t]=link.x[No]
        robot.storage.q[ind][t]=link.q[No]
    end
    for (ind,constraint) in enumerate(robot.constraints)
        robot.storage.λ[ind][t]=constraint.s1
    end
end

@inline function updatePos!(link::Link,dt)
    x2 = link.x[2]
    q2 = link.q[2]
    link.x[1] = x2
    link.x[2] = x2 + getvnew(link)*dt
    link.q[1] = q2
    link.q[2] = dt/2*(Lmat(q2)*ωbar(link,dt))
    return
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
            θ[i,t] = angleaxis(q)[1]*sign(angleaxis(q)[2][1])
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
