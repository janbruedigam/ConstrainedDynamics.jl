mutable struct Mechanism{T,N}
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
    function Mechanism(origin::Origin{T},links::Vector{Link{T}},constraints::Vector{<:Constraint{T}};
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

    function Mechanism(origin::Origin{T},links::Vector{Link{T}};
        tend::T=10., dt::T=.01, g::T=-9.81, No=2) where T

        constraints = Vector{Constraint{T}}(undef,0)
        for link in links
            push!(constraints,Constraint(OriginConnection(origin,link)))
        end
        Mechanism(origin,links,constraints,tend=tend, dt=dt, g=g, No=No)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, R::Mechanism{T}) where {T}
    summary(io, R); println(io, " with ", length(R.links), " links and ", length(R.constraints), " constraints")
end

function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id,link) in pairs(mechanism.links)
        for cid in directchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)),id,getconstraint(mechanism,cid),mechanism)
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,link,mechanism)
    end

    for node in mechanism.constraints
        id = node.id

        for cid in directchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)),node,cid,mechanism)
        end

        for cid in loopchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)))
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,node,mechanism)
    end
end

@inline getlink(mechanism::Mechanism,id::Int64) = mechanism.links[id]
@inline getlink(mechanism::Mechanism,id::Nothing) = mechanism.origin
@inline getconstraint(mechanism::Mechanism,id::Int64) = mechanism.constraints[id]

# @inline function getnode(mechanism::Mechanism,id::Int64) # should only be used in setup
#      if haskey(mechanism.ldict,id)
#          return getlink(mechanism,id)
#      elseif haskey(mechanism.cdict,id)
#          return getconstraint(mechanism,id)
#      elseif id == mechanism.originid
#          return mechanism.origin
#      else
#          error("not found.")
#      end
#  end

@inline function normf(link::Link{T},mechanism::Mechanism) where T
    f = dynamics(link,mechanism)
    return dot(f,f)
end

@inline function normf(c::Constraint,mechanism::Mechanism)
    f = g(c,mechanism)
    return dot(f,f)
end

@inline function GtλTof!(link::Link,c::Constraint,mechanism)
    link.f -= ∂g∂pos(c,link.id,mechanism)'*c.s1
    return
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0

    for link in mechanism.links
        mechanism.normf += normf(link,mechanism)
    end
    foreach(addNormf!,mechanism.constraints,mechanism)

    return sqrt(mechanism.normf)
end

@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    mechanism.normΔs += mapreduce(normΔs,+,mechanism.links)
    foreach(addNormΔs!,mechanism.constraints,mechanism)

    return sqrt(mechanism.normΔs)
end

@inline function addNormf!(c::Constraint,mechanism::Mechanism)
    mechanism.normf += normf(c,mechanism)
    return
end

@inline function addNormΔs!(component::Component,mechanism::Mechanism)
    mechanism.normΔs += normΔs(component)
    return
end

function saveToTraj!(mechanism::Mechanism,t)
    No = mechanism.No
    for (ind,link) in enumerate(mechanism.links)
        mechanism.storage.x[ind][t]=link.x[No]
        mechanism.storage.q[ind][t]=link.q[No]
    end
    for (ind,constraint) in enumerate(mechanism.constraints)
        mechanism.storage.λ[ind][t]=constraint.s1
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


function simulate!(mechanism::Mechanism;save::Bool=false,debug::Bool=false,disp::Bool=false)
    links = mechanism.links
    constraints = mechanism.constraints
    dt = mechanism.dt
    foreach(s0tos1!,links)
    foreach(s0tos1!,constraints)

    for i=mechanism.steps
        newton!(mechanism,warning=debug)
        save && saveToTraj!(mechanism,i)
        foreach(updatePos!,links,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return
end



function plotθ(mechanism::Mechanism{T},id) where T
    n = length(mechanism.links)
    θ = zeros(T,n,length(mechanism.steps))
    for i=1:n
        qs = mechanism.storage.q[i]
        for (t,q) in enumerate(qs)
            θ[i,t] = angleaxis(q)[1]*sign(angleaxis(q)[2][1])
        end
    end

    p = plot(collect(0:mechanism.dt:mechanism.tend-mechanism.dt),θ[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:mechanism.dt:mechanism.tend-mechanism.dt),θ[ind,:])
    end
    return p
end

function plotλ(mechanism::Mechanism{T},id) where T
    n = sum(length.(mechanism.constraints))
    λ = zeros(T,n,length(mechanism.steps))
    startpos = 1
    endpos = 0
    for i=1:length(mechanism.constraints)
        endpos = startpos + length(mechanism.constraints[i]) -1

        λs = mechanism.storage.λ[i]
        for (t,val) in enumerate(λs)
            λ[startpos:endpos,t] = val
        end

        startpos = endpos + 1
    end

    p = plot(collect(0:mechanism.dt:mechanism.tend-mechanism.dt),λ[id[1],:])
    for ind in Iterators.rest(id,2)
        plot!(collect(0:mechanism.dt:mechanism.tend-mechanism.dt),λ[ind,:])
    end
    return p
end
