mutable struct Mechanism{T,N}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T
    No::Int64

    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    #TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T

    graph::Graph{N}

    ldu::SparseLDU{T}
    storage::Storage{T}

    μ::Float64
    αmax::Float64

    #TODO no constraints input
    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},
        eqconstraints::Vector{<:EqualityConstraint{T}}, ineqconstraints::Vector{<:InequalityConstraint{T}};
        tend::T=10., dt::T=.01, g::T=-9.81, No=2) where T


        resetGlobalID()

        Nb = length(bodies)
        Ne = length(eqconstraints)
        Ni = length(ineqconstraints)
        N = Nb+Ne
        steps = Int(ceil(tend/dt))

        currentid = 1

        bdict = Dict{Int64,Int64}()
        for (ind,body) in enumerate(bodies)
            push!(body.x, [body.x[1] for i=1:No-1]...)
            push!(body.q, [body.q[1] for i=1:No-1]...)
            push!(body.F, [body.F[1] for i=1:No-1]...)
            push!(body.τ, [body.τ[1] for i=1:No-1]...)

            for c in eqconstraints
                c.pid == body.id && (c.pid = currentid)
                for (ind,bodyid) in enumerate(c.bodyids)
                    if bodyid == body.id
                        c.bodyids = setindex(c.bodyids,currentid,ind)
                        c.constraints[ind].cid = currentid
                    end
                end
            end

            for c in ineqconstraints
                c.pid == body.id && (c.pid = currentid)
            end

            body.id = currentid
            currentid+=1

            bdict[body.id] = ind
        end

        eqdict = Dict{Int64,Int64}()
        for (ind,c) in enumerate(eqconstraints)
            c.id = currentid
            currentid+=1

            eqdict[c.id] = ind
        end

        ineqdict = Dict{Int64,Int64}()
        for (ind,c) in enumerate(ineqconstraints)
            c.id = currentid
            currentid+=1

            ineqdict[c.id] = ind
        end

        normf = zero(T)
        normΔs = zero(T)

        graph = Graph(origin,bodies,eqconstraints,ineqconstraints)
        ldu = SparseLDU(graph,bodies,eqconstraints,ineqconstraints,bdict,eqdict,ineqdict)

        storage = Storage{T}(steps,Nb,Ne)

        bodies = UnitDict(bodies)
        eqconstraints = UnitDict((eqconstraints[1].id):(eqconstraints[Ne].id),eqconstraints)
        if !isempty(ineqconstraints)
            ineqconstraints = UnitDict((ineqconstraints[1].id):(ineqconstraints[Ni].id),ineqconstraints)
        else
            ineqconstraints = UnitDict(0:0,ineqconstraints)
        end
        new{T,N}(tend,Base.OneTo(steps),dt,g,No,origin,bodies,eqconstraints,ineqconstraints,normf,normΔs,graph,ldu,storage,1,1)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}};
        tend::T=10., dt::T=.01, g::T=-9.81, No=2) where T

        constraints = Vector{EqualityConstraint{T}}(undef,0)
        for body in bodies
            push!(constraints,EqualityConstraint(OriginConnection(origin,body)))
        end
        Mechanism(origin,bodies,constraints,tend=tend, dt=dt, g=g, No=No)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},constraints::Vector{<:EqualityConstraint{T}};
        tend::T=10., dt::T=.01, g::T=-9.81, No=2) where T

        ineqconstraints = Vector{InequalityConstraint{T}}(undef,0)
        Mechanism(origin,bodies,constraints,ineqconstraints,tend=tend, dt=dt, g=g, No=No)
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, M::Mechanism{T}) where {T}
    summary(io, M); println(io, " with ", length(M.bodies), " bodies and ", length(M.eqconstraints), " constraints")
end

function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id,body) in pairs(mechanism.bodies)
        for cid in directchildren(graph,id)
            setLU!(getentry(ldu,(id,cid)),id,geteqconstraint(mechanism,cid),mechanism)
        end

        diagonal = getentry(ldu,id)
        setDandŝ!(diagonal,body,mechanism)
        for cid in ineqchildren(graph,id)
            extendDandŝ!(diagonal,body,getineqconstraint(mechanism,cid),mechanism)
        end
    end

    for node in mechanism.eqconstraints
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

@inline getbody(mechanism::Mechanism,id::Int64) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism,id::Nothing) = mechanism.origin
@inline geteqconstraint(mechanism::Mechanism,id::Int64) = mechanism.eqconstraints[id]
@inline getineqconstraint(mechanism::Mechanism,id::Int64) = mechanism.ineqconstraints[id]

# @inline function getnode(mechanism::Mechanism,id::Int64) # should only be used in setup
#      if haskey(mechanism.bdict,id)
#          return getbody(mechanism,id)
#      elseif haskey(mechanism.cdict,id)
#          return getconstraint(mechanism,id)
#      elseif id == mechanism.originid
#          return mechanism.origin
#      else
#          error("not found.")
#      end
#  end

@inline function normf(body::Body{T},mechanism::Mechanism) where T
    f = dynamics(body,mechanism)
    return dot(f,f)
end

@inline function normf(c::EqualityConstraint,mechanism::Mechanism)
    f = g(c,mechanism)
    return dot(f,f)
end

@inline function normf(c::InequalityConstraint,mechanism::Mechanism)
    f = g(c,mechanism)
    d = h(c,mechanism)
    return dot([f;d],[f;d])
end

@inline function normfμ(c::InequalityConstraint,mechanism::Mechanism)
    f = g(c,mechanism)
    d = hμ(c,mechanism)
    return dot([f;d],[f;d])
end

@inline function GtλTof!(body::Body,c::EqualityConstraint,mechanism)
    body.f -= ∂g∂pos(c,body.id,mechanism)'*c.s1
    return
end

@inline function NtγTof!(body::Body,c::InequalityConstraint,mechanism)
    Nx = SVector{6,Float64}(0,0,1,0,0,0)'
    body.f -= Nx'*c.ga1

    D = [1 0 0 0 0 0;0 1 0 0 0 0]
    body.f -= D'*c.b1
    return
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0

    for body in mechanism.bodies
        mechanism.normf += normf(body,mechanism)
    end
    foreach(addNormf!,mechanism.eqconstraints,mechanism)
    for ineq in mechanism.ineqconstraints
        mechanism.normf += normf(ineq,mechanism)
    end

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    mechanism.normf = 0

    for body in mechanism.bodies
        mechanism.normf += normf(body,mechanism)
    end
    foreach(addNormf!,mechanism.eqconstraints,mechanism)
    # mechanism.normf = sqrt(mechanism.normf)
    for ineq in mechanism.ineqconstraints
        mechanism.normf += normfμ(ineq,mechanism)
        # mechanism.normf -= mechanism.μ*ineq.sl1
    end

    # return sqrt(mechanism.normf)
    return mechanism.normf
end

@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    mechanism.normΔs += mapreduce(normΔs,+,mechanism.bodies)
    foreach(addNormΔs!,mechanism.eqconstraints,mechanism)
    foreach(addNormΔs!,mechanism.ineqconstraints,mechanism)

    return sqrt(mechanism.normΔs)
end

@inline function addNormf!(c::EqualityConstraint,mechanism::Mechanism)
    mechanism.normf += normf(c,mechanism)
    return
end

@inline function addNormΔs!(component::Component,mechanism::Mechanism)
    mechanism.normΔs += normΔs(component)
    return
end

function computeα!(mechanism::Mechanism)
    ldu = mechanism.ldu

    τ = 0.995
    αmax = 1.

    for ineq in mechanism.ineqconstraints
        sl = getineq(ldu,ineq.id).sl
        ga = getineq(ldu,ineq.id).ga
        slf = getineq(ldu,ineq.id).slf
        psi = getineq(ldu,ineq.id).psi

        if sl > 0
            temp = minimum([1.;τ*ineq.sl1/sl])
            αmax = minimum([αmax;temp])
        end

        if ga > 0
            temp = minimum([1.;τ*ineq.ga1/ga])
            αmax = minimum([αmax;temp])
        end

        if slf > 0
            temp = minimum([1.;τ*ineq.slf1/slf])
            αmax = minimum([αmax;temp])
        end

        if psi > 0
            temp = minimum([1.;τ*ineq.psi1/psi])
            αmax = minimum([αmax;temp])
        end
    end

    mechanism.αmax = αmax

    return
end

function saveToTraj!(mechanism::Mechanism,t)
    No = mechanism.No
    for (ind,body) in enumerate(mechanism.bodies)
        mechanism.storage.x[ind][t]=body.x[No]
        mechanism.storage.q[ind][t]=body.q[No]
    end
    for (ind,constraint) in enumerate(mechanism.eqconstraints)
        mechanism.storage.λ[ind][t]=constraint.s1
    end
end

@inline function updatePos!(body::Body,dt)
    x2 = body.x[2]
    q2 = body.q[2]
    body.x[1] = x2
    body.x[2] = x2 + getvnew(body)*dt
    body.q[1] = q2
    body.q[2] = dt/2*(Lmat(q2)*ωbar(body,dt))
    return
end


function simulate!(mechanism::Mechanism;save::Bool=false,debug::Bool=false,disp::Bool=false)
    bodies = mechanism.bodies
    constraints = mechanism.eqconstraints
    dt = mechanism.dt
    foreach(s0tos1!,bodies)
    foreach(s0tos1!,constraints)

    for i=mechanism.steps
        newton!(mechanism,warning=debug)
        save && saveToTraj!(mechanism,i)
        foreach(updatePos!,bodies,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return
end

function simulate_ip!(mechanism::Mechanism;save::Bool=false,debug::Bool=false,disp::Bool=false)
    bodies = mechanism.bodies
    eqconstraints = mechanism.eqconstraints
    ineqconstraints = mechanism.ineqconstraints
    dt = mechanism.dt
    foreach(s0tos1!,bodies)
    foreach(s0tos1!,eqconstraints)
    foreach(s0tos1!,ineqconstraints)

    for i=mechanism.steps
        # newton!(mechanism,warning=debug)
        # newton_ip!(mechanism,bodies[1])
        newton_ip!(mechanism,warning=debug)
        save && saveToTraj!(mechanism,i)
        foreach(updatePos!,bodies,dt)

        disp && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return
end


function plotθ(mechanism::Mechanism{T},id) where T
    n = length(mechanism.bodies)
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
    n = sum(length.(mechanism.eqconstraints))
    λ = zeros(T,n,length(mechanism.steps))
    startpos = 1
    endpos = 0
    for i=1:length(mechanism.eqconstraints)
        endpos = startpos + length(mechanism.eqconstraints[i]) -1

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
