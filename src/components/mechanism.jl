mutable struct Mechanism{T,N,Ni}
    tend::T
    steps::Base.OneTo{Int64}
    Δt::T
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

    shapes::Vector{<:Shape{T}}


    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},
        eqcs::Vector{<:EqualityConstraint{T}}, ineqcs::Vector{<:InequalityConstraint{T}};
        tend::T = 10., Δt::T = .01, g::T = -9.81, No = 2, shapes::Vector{<:Shape{T}} = Shape{T}[]) where T


        resetGlobalID()

        Nb = length(bodies)
        Ne = length(eqcs)
        Ni = length(ineqcs)
        N = Nb + Ne
        steps = Int(ceil(tend / Δt))

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

            for shape in shapes
                for (i,id) in enumerate(shape.bodyids)
                    id == body.id && (shape.bodyids[i] = currentid)
                end
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

        new{T,N,Ni}(tend, Base.OneTo(steps), Δt, g, No, origin, bodies, eqcs, ineqcs, normf, normΔs, graph, ldu, storage, α, μ, shapes)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},eqcs::Vector{<:EqualityConstraint{T}};
        tend::T = 10., Δt::T = .01, g::T = -9.81, No = 2, shapes::Vector{<:Shape{T}} = Shape{T}[]) where T

        ineqcs = InequalityConstraint{T}[]
        Mechanism(origin, bodies, eqcs, ineqcs, tend = tend, Δt = Δt, g = g, No = No, shapes = shapes)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},ineqcs::Vector{<:InequalityConstraint{T}};
        tend::T = 10., Δt::T = .01, g::T = -9.81, No = 2, shapes::Vector{<:Shape{T}} = Shape{T}[]) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(OriginConnection(origin, body)))
        end
        Mechanism(origin, bodies, eqc, ineqcs, tend = tend, Δt = Δt, g = g, No = No, shapes = shapes)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}};
        tend::T = 10., Δt::T = .01, g::T = -9.81, No = 2, shapes::Vector{<:Shape{T}} = Shape{T}[]) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(OriginConnection(origin, body)))
        end
        Mechanism(origin, bodies, eqc, tend = tend, Δt = Δt, g = g, No = No, shapes = shapes)
    end    

    function Mechanism(filename::AbstractString; scalar_type::Type{T}=Float64, tend::T=10., Δt::T=.01, g::T=-9.81, No::Int64=2) where T
        xdoc = LightXML.parse_file(filename)
        xroot = LightXML.root(xdoc)
        @assert LightXML.name(xroot) == "robot"
    
        xlinks = get_elements_by_tagname(xroot, "link")
        xjoints = get_elements_by_tagname(xroot, "joint")
    
        ldict = parse_links(xlinks,T)
    
        origin, links, joints = parse_joints(xjoints,ldict)
    
        free(xdoc)
    
        Mechanism(origin, links, joints, tend=tend, Δt=Δt, g=g, No=No)
        # return origin, links, joints
        
    end
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, M::Mechanism{T,N,0}) where {T,N}
    summary(io, M); println(io, " with ", length(M.bodies), " bodies and ", length(M.eqconstraints), " constraints")
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, M::Mechanism{T,N,Ni}) where {T,N,Ni}
    summary(io, M); println(io, " with ", length(M.bodies), " bodies, ", length(M.eqconstraints), " equality constraints, and ", Ni, " inequality constraints")
end


@inline getbody(mechanism::Mechanism, id::Int64) = mechanism.bodies[id]
@inline getbody(mechanism::Mechanism, id::Nothing) = mechanism.origin
@inline geteqconstraint(mechanism::Mechanism, id::Int64) = mechanism.eqconstraints[id]
@inline getineqconstraint(mechanism::Mechanism, id::Int64) = mechanism.ineqconstraints[id]


function setPosition!(mechanism::Mechanism{T}, body::Body{T};x::AbstractVector{T}=SVector{3,T}(0,0,0),q::Quaternion{T}=Quaternion{T}()) where T
    for i=1:mechanism.No
        body.x[i] = x
        body.q[i] = q
    end
end

function setPosition!(mechanism::Mechanism{T}, body1::Body{T}, body2::Body{T};
    p1::AbstractVector{T}=SVector{3,T}(0,0,0), p2::AbstractVector{T}=SVector{3,T}(0,0,0), Δx::AbstractVector{T}=SVector{3,T}(0,0,0),Δq::Quaternion{T}=Quaternion{T}()) where T

    q = body1.q[1]*Δq
    x = body1.x[1] + vrotate(SVector{3,T}(p1+Δx), body1.q[1]) - vrotate(SVector{3,T}(p2), q)
    
    setPosition!(mechanism, body2;x=x,q=q)
end

function setPosition!(mechanism::Mechanism{T}, body1::Origin{T}, body2::Body{T};
    p1::AbstractVector{T}=SVector{3,T}(0,0,0), p2::AbstractVector{T}=SVector{3,T}(0,0,0), Δx::AbstractVector{T}=SVector{3,T}(0,0,0),Δq::Quaternion{T}=Quaternion{T}()) where T

    q = Δq
    x = p1+Δx - vrotate(SVector{3,T}(p2), q)
    

    setPosition!(mechanism, body2;x=x,q=q)
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body::Body{T};v::AbstractVector{T}=SVector{3,T}(0,0,0),ω::AbstractVector{T}=SVector{3,T}(0,0,0)) where T
    body.s0 = [v;ω]
    s0tos1!(body)
    updatePos!(body, mechanism.Δt)
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body1::Body{T}, body2::Body{T};
    p1::AbstractVector{T}=SVector{3,T}(0,0,0), p2::AbstractVector{T}=SVector{3,T}(0,0,0), Δv2::AbstractVector{T}=SVector{3,T}(0,0,0),Δω2::AbstractVector{T}=SVector{3,T}(0,0,0)) where T

    Δt = mechanism.Δt
    x2 = body1.x[2] + vrotate(SVector{3,T}(p1), body1.q[2]) - vrotate(SVector{3,T}(p2), body2.q[2]) + vrotate(SVector{3,T}(Δv2),body1.q[1])*Δt
    q2 = body1.q[2]*(Δt/2 * Quaternion(sqrt(4 / Δt^2 - dot(Δω2, Δω2)), SVector{3,T}(Δω2)))

    v = (x2-body2.x[1])/Δt
    ω = 2/Δt*Vmat(body2.q[1]\q2)

    setVelocity!(mechanism, body2;v=v,ω=ω)
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body1::Origin{T}, body2::Body{T};
    p1::AbstractVector{T}=SVector{3,T}(0,0,0), p2::AbstractVector{T}=SVector{3,T}(0,0,0), Δv2::AbstractVector{T}=SVector{3,T}(0,0,0),Δω2::AbstractVector{T}=SVector{3,T}(0,0,0)) where T

    Δt = mechanism.Δt
    x2 = p1 - vrotate(SVector{3,T}(p2), body2.q[2]) + Δv2*Δt
    q2 = Δt/2 * Quaternion(sqrt(4 / Δt^2 - dot(Δω2, Δω2)), SVector{3,T}(Δω2))

    v = (x2-body2.x[1])/Δt
    ω = 2/Δt*Vmat(body2.q[1]\q2)

    setVelocity!(mechanism, body2;v=v,ω=ω)
end

function setForce!(mechanism::Mechanism{T}, body::Body{T};F::AbstractVector{T}=SVector{3,T}(0,0,0),r::AbstractVector{T}=SVector{3,T}(0,0,0),τ::AbstractVector{T}=SVector{3,T}(0,0,0)) where T
    τ += torqueFromForce(F,r)
    setForce!(body,F,τ,mechanism.No)
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

    p = plot(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), θ[id[1],:])
    for ind in Iterators.rest(id, 2)
        plot!(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), θ[ind,:])
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

    p = plot(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), λ[id[1],:])
    for ind in Iterators.rest(id, 2)
        plot!(collect(0:mechanism.Δt:mechanism.tend - mechanism.Δt), λ[ind,:])
    end
    return p
end
