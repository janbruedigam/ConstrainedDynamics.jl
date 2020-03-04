mutable struct Mechanism{T,N,Ni}
    tend::T
    steps::Base.OneTo{Int64}
    dt::T
    g::T
    No::Int64 # order of integrator, currently only No=2 (1st order) implemented

    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    # TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T

    graph::Graph{N}

    ldu::SparseLDU{T}
    storage::Storage{T}

    α::T
    μ::T


    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},
        eqcs::Vector{<:EqualityConstraint{T}}, ineqcs::Vector{<:InequalityConstraint{T}};
        tend::T = 10., dt::T = .01, g::T = -9.81, No = 2) where T


        resetGlobalID()

        Nb = length(bodies)
        Ne = length(eqcs)
        Ni = length(ineqcs)
        N = Nb + Ne
        steps = Int(ceil(tend / dt))

        currentid = 1

        bdict = Dict{Int64,Int64}()
        for (ind, body) in enumerate(bodies)
            push!(body.x, [body.x[1] for i = 1:No - 1]...)
            push!(body.q, [body.q[1] for i = 1:No - 1]...)
            push!(body.F, [body.F[1] for i = 1:No - 1]...)
            push!(body.τ, [body.τ[1] for i = 1:No - 1]...)

            for eqc in eqcs
                eqc.pid == body.id && (eqc.pid = currentid)
                for (ind, bodyid) in enumerate(eqc.bodyids)
                    if bodyid == body.id
                        eqc.bodyids = setindex(eqc.bodyids, currentid, ind)
                        eqc.constraints[ind].cid = currentid
                    end
                end
            end

            for ineqc in ineqcs
                ineqc.pid == body.id && (ineqc.pid = currentid)
            end

            body.id = currentid
            currentid += 1

            bdict[body.id] = ind
        end

        eqdict = Dict{Int64,Int64}()
        for (ind, eqc) in enumerate(eqcs)
            eqc.id = currentid
            currentid += 1

            eqdict[eqc.id] = ind
        end

        ineqdict = Dict{Int64,Int64}()
        for (ind, ineqc) in enumerate(ineqcs)
            ineqc.id = currentid
            currentid += 1

            ineqdict[ineqc.id] = ind
        end

        normf = 0
        normΔs = 0

        graph = Graph(origin, bodies, eqcs, ineqcs)
        ldu = SparseLDU(graph, bodies, eqcs, ineqcs, bdict, eqdict, ineqdict)

        storage = Storage{T}(steps, Nb, Ne)

        bodies = UnitDict(bodies)
        eqcs = UnitDict((eqcs[1].id):(eqcs[Ne].id), eqcs)
        if Ni > 0
            ineqcs = UnitDict((ineqcs[1].id):(ineqcs[Ni].id), ineqcs)
        else
            ineqcs = UnitDict(0:0, ineqcs)
        end

        α = 1
        μ = 1

        new{T,N,Ni}(tend, Base.OneTo(steps), dt, g, No, origin, bodies, eqcs, ineqcs, normf, normΔs, graph, ldu, storage, α, μ)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},eqcs::Vector{<:EqualityConstraint{T}};
        tend::T = 10., dt::T = .01, g::T = -9.81, No = 2) where T

        ineqcs = InequalityConstraint{T}[]
        Mechanism(origin, bodies, eqcs, ineqcs, tend = tend, dt = dt, g = g, No = No)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},ineqcs::Vector{<:InequalityConstraint{T}};
        tend::T = 10., dt::T = .01, g::T = -9.81, No = 2) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(OriginConnection(origin, body)))
        end
        Mechanism(origin, bodies, eqc, ineqcs, tend = tend, dt = dt, g = g, No = No)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}};
        tend::T = 10., dt::T = .01, g::T = -9.81, No = 2) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(OriginConnection(origin, body)))
        end
        Mechanism(origin, bodies, eqc, tend = tend, dt = dt, g = g, No = No)
    end    
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, M::Mechanism{T,N,0}) where {T,N}
    summary(io, M); println(io, " with ", length(M.bodies), " bodies and ", length(M.eqconstraints), " constraints")
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, M::Mechanism{T,N,Ni}) where {T,N,Ni}
    summary(io, M); println(io, " with ", length(M.bodies), " bodies, ", length(M.eqconstraints), " equality constraints, and ", Ni, " inequality constraints")
end

function setentries!(mechanism::Mechanism)
    graph = mechanism.graph
    ldu = mechanism.ldu

    for (id, body) in pairs(mechanism.bodies)
        for cid in directchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)), id, geteqconstraint(mechanism, cid), mechanism)
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(diagonal, body, mechanism)
        for cid in ineqchildren(graph, id)
            extendDandΔs!(diagonal, body, getineqconstraint(mechanism, cid), mechanism)
        end
    end

    for node in mechanism.eqconstraints
        id = node.id

        for cid in directchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)), node, cid, mechanism)
        end

        for cid in loopchildren(graph, id)
            setLU!(getentry(ldu, (id, cid)))
        end

        diagonal = getentry(ldu, id)
        setDandΔs!(diagonal, node, mechanism)
    end
end

@inline getbody(mechanism::Mechanism, id::Int64) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism, id::Nothing) = mechanism.origin
@inline geteqconstraint(mechanism::Mechanism, id::Int64) = mechanism.eqconstraints[id]
@inline getineqconstraint(mechanism::Mechanism, id::Int64) = mechanism.ineqconstraints[id]


@inline function normf(body::Body{T}, mechanism::Mechanism) where T
    f = dynamics(body, mechanism)
    return dot(f, f)
end

@inline function normf(c::EqualityConstraint, mechanism::Mechanism)
    f = g(c, mechanism)
    return dot(f, f)
end

@inline function normf(ineqc::InequalityConstraint, mechanism::Mechanism)
    f = gs(ineqc, mechanism)
    d = h(ineqc)
    return dot(f, f) + dot(d, d)
end

@inline function normfμ(ineqc::InequalityConstraint, mechanism::Mechanism)
    f = gs(ineqc, mechanism)
    d = hμ(ineqc, mechanism.μ)
    return dot(f, f) + dot(d, d)
end

@inline function GtλTof!(body::Body, eqc::EqualityConstraint, mechanism)
    body.f -= ∂g∂pos(eqc, body.id, mechanism)' * eqc.s1
    return
end

@inline function NtγTof!(body::Body, ineqc::InequalityConstraint, mechanism)
    body.f -= ∂g∂pos(ineqc, body, mechanism)' * ineqc.γ1
    return
end

@inline function normf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(body, mechanism)
    end
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormf!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function meritf(mechanism::Mechanism)
    mechanism.normf = 0

    # Allocates otherwise
    for body in mechanism.bodies
        mechanism.normf += normf(body, mechanism)
    end
    foreach(addNormf!, mechanism.eqconstraints, mechanism)
    foreach(addNormfμ!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normf)
end

@inline function normΔs(mechanism::Mechanism)
    mechanism.normΔs = 0

    # Allocates otherwise
    mechanism.normΔs += mapreduce(normΔs, +, mechanism.bodies)
    foreach(addNormΔs!, mechanism.eqconstraints, mechanism)
    foreach(addNormΔs!, mechanism.ineqconstraints, mechanism)

    return sqrt(mechanism.normΔs)
end

@inline function addNormf!(ineqc::InequalityConstraint, mechanism::Mechanism)
    mechanism.normf += normf(ineqc, mechanism)
    return
end

@inline function addNormfμ!(ineqc::InequalityConstraint, mechanism::Mechanism)
    mechanism.normf += normfμ(ineqc, mechanism)
    return
end

@inline function addNormf!(eqc::EqualityConstraint, mechanism::Mechanism)
    mechanism.normf += normf(eqc, mechanism)
    return
end

@inline function addNormΔs!(component::Component, mechanism::Mechanism)
    mechanism.normΔs += normΔs(component)
    return
end

function feasibilityStepLength!(mechanism::Mechanism)
    ldu = mechanism.ldu

    τ = 0.995
    mechanism.α = 1.

    for ineqc in mechanism.ineqconstraints
        feasibilityStepLength!(ineqc, getineqentry(ldu, ineqc.id), τ, mechanism)
    end

    return
end

function feasibilityStepLength!(ineqc::InequalityConstraint{T,N}, ineqentry::InequalityEntry, τ, mechanism) where {T,N}
    s1 = ineqc.s1
    γ1 = ineqc.γ1
    Δs = ineqentry.Δs
    Δγ = ineqentry.Δγ

    for i = 1:N
        αmax = τ * s1[i] / Δs[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
        αmax = τ * γ1[i] / Δγ[i]
        (αmax > 0) && (αmax < mechanism.α) && (mechanism.α = αmax)
    end
    return
end


function saveToStorage!(mechanism::Mechanism, t)
    No = mechanism.No
    for (ind, body) in enumerate(mechanism.bodies)
        mechanism.storage.x[ind][t] = body.x[No]
        mechanism.storage.q[ind][t] = body.q[No]
    end
    for (ind, eqc) in enumerate(mechanism.eqconstraints)
        mechanism.storage.eqmultipliers[ind][t] = eqc.s1
    end
end

@inline function updatePos!(body::Body, dt)
    body.x[1] = body.x[2]
    body.x[2] = getx3(body, dt)
    body.q[1] = body.q[2]
    body.q[2] = getq3(body, dt)
    return
end


function simulate!(mechanism::Mechanism;save::Bool = false,debug::Bool = false)
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints
    ineqcs = mechanism.ineqconstraints
    dt = mechanism.dt
    foreach(s0tos1!, bodies)
    foreach(s0tos1!, eqcs)
    foreach(s0tos1!, ineqcs)

    for i = mechanism.steps
        newton!(mechanism, warning = debug)
        save && saveToStorage!(mechanism, i)
        foreach(updatePos!, bodies, dt)

        # debug && (i*dt)%1<dt*(1.0-.1) && display(i*dt)
    end
    return
end


function plotθ(mechanism::Mechanism{T}, id) where T
    n = length(mechanism.bodies)
    θ = zeros(T, n, length(mechanism.steps))
    for i = 1:n
        qs = mechanism.storage.q[i]
        for (t, q) in enumerate(qs)
            θ[i,t] = angleaxis(q)[1] * sign(angleaxis(q)[2][1])
        end
    end

    p = plot(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), θ[id[1],:])
    for ind in Iterators.rest(id, 2)
        plot!(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), θ[ind,:])
    end
    return p
end

function plotλ(mechanism::Mechanism{T}, id) where T
    n = sum(length.(mechanism.eqconstraints))
    λ = zeros(T, n, length(mechanism.steps))
    startpos = 1
    endpos = 0
    for i = 1:length(mechanism.eqconstraints)
        endpos = startpos + length(mechanism.eqconstraints[i]) - 1

        λs = mechanism.storage.λ[i]
        for (t, val) in enumerate(λs)
            λ[startpos:endpos,t] = val
        end

        startpos = endpos + 1
    end

    p = plot(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), λ[id[1],:])
    for ind in Iterators.rest(id, 2)
        plot!(collect(0:mechanism.dt:mechanism.tend - mechanism.dt), λ[ind,:])
    end
    return p
end
