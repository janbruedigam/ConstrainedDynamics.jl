abstract type AbstractMechanism{T,Nn,Nb,Ne,Ni} end

mutable struct Mechanism{T,Nn,Nb,Ne,Ni} <: AbstractMechanism{T,Nn,Nb,Ne,Ni}
    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    graph::Graph{Nn}
    ldu::SparseLDU{T}

    # TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T

    Δt::T
    g::T

    α::T
    μ::T


    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},
            eqcs::Vector{<:EqualityConstraint{T}}, ineqcs::Vector{<:InequalityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81, shapes::Vector{<:Shape{T}} = Shape{T}[]
        ) where T

        resetGlobalID()
        order = getGlobalOrder()

        for body in bodies
            if norm(body.m) == 0 || norm(body.J) == 0
                @info "Potentially bad inertial properties detected"
            end 
        end

        Nb = length(bodies)
        Ne = length(eqcs)
        Ni = length(ineqcs)
        Nn = Nb + Ne
        # steps = Int64(ceil(tend / Δt))

        currentid = 1

        bdict = Dict{Int64,Int64}()
        for (ind, body) in enumerate(bodies)
            initknotpoints!(body.state, order)

            for eqc in eqcs
                eqc.parentid == body.id && (eqc.parentid = currentid)
                for (ind, bodyid) in enumerate(eqc.childids)
                    if bodyid == body.id
                        eqc.childids = setindex(eqc.childids, currentid, ind)
                        eqc.constraints[ind].childid = currentid
                    end
                end
            end

            for ineqc in ineqcs
                ineqc.parentid == body.id && (ineqc.parentid = currentid)
            end

            for shape in shapes
                for (i, id) in enumerate(shape.bodyids)
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

        # storage = Storage{T}(steps, Nb, Ne)

        bodies = UnitDict(bodies)
        eqcs = UnitDict((eqcs[1].id):(eqcs[Ne].id), eqcs)
        if Ni > 0
            ineqcs = UnitDict((ineqcs[1].id):(ineqcs[Ni].id), ineqcs)
        else
            ineqcs = UnitDict(0:0, ineqcs)
        end

        α = 1
        μ = 1

        new{T,Nn,Nb,Ne,Ni}(origin, bodies, eqcs, ineqcs, graph, ldu, normf, normΔs, Δt, g, α, μ)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},eqcs::Vector{<:EqualityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81, shapes::Vector{<:Shape{T}} = Shape{T}[]
        ) where T

        ineqcs = InequalityConstraint{T}[]
        return Mechanism(origin, bodies, eqcs, ineqcs, Δt = Δt, g = g, shapes = shapes)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},ineqcs::Vector{<:InequalityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81, shapes::Vector{<:Shape{T}} = Shape{T}[]
        ) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(Floating(origin, body)))
        end
        return Mechanism(origin, bodies, eqc, ineqcs, Δt = Δt, g = g, shapes = shapes)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}};
            Δt::Real = .01, g::Real = -9.81, shapes::Vector{<:Shape{T}} = Shape{T}[]
        ) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(Floating(origin, body)))
        end
        return Mechanism(origin, bodies, eqc, Δt = Δt, g = g, shapes = shapes)
    end

    function Mechanism(filename::AbstractString; floating::Bool=false, type::Type{T} = Float64, Δt::Real = .01, g::Real = -9.81) where T
        origin, links, joints, shapes = parse_urdf(filename, floating, T)

        mechanism = Mechanism(origin, links, joints, shapes = shapes, Δt = Δt, g = g)

        graph = mechanism.graph
        xjointlist = Dict{Int64,SVector{3,T}}() # stores id, x in world frame
        qjointlist = Dict{Int64,UnitQuaternion{T}}() # stores id, q in world frame

        for id in graph.rdfslist # from root to leaves
            component = getcomponent(mechanism, id)
            if component isa Body
                body = component
                xbodylocal = body.state.xc
                qbodylocal = body.state.qc
                shape = getshape(shapes, id)

                preds = predecessors(graph, id);                    @assert length(preds) == 1
                parentid = preds[1]
                constraint = geteqconstraint(mechanism, parentid);  @assert length(constraint.constraints) == 2
                grandpreds = predecessors(graph, parentid);         @assert length(grandpreds) ∈ [0;1]

                if length(grandpreds) == 1 # predecessor is link
                    grandparentid = grandpreds[1]

                    parentbody = getbody(mechanism, grandparentid)
                    grandgrandpreds = predecessors(graph, grandparentid);               @assert length(grandgrandpreds) == 1
                    grandgrandparentid = grandgrandpreds[1]
                    parentconstraint = geteqconstraint(mechanism, grandgrandparentid);  @assert length(parentconstraint.constraints) == 2

                    xparentbody = parentbody.state.xc # in world frame
                    qparentbody = parentbody.state.qc # in world frame

                    xparentjoint = xjointlist[parentconstraint.id] # in world frame
                    qparentjoint = qjointlist[parentconstraint.id] # in world frame
                else # predecessor is origin
                    parentbody = origin

                    xparentbody = SA{T}[0; 0; 0]
                    qparentbody = one(UnitQuaternion{T})

                    xparentjoint = SA{T}[0; 0; 0]
                    qparentjoint = one(UnitQuaternion{T})
                end

                # urdf joint's x and q in parent's (parentbody) frame
                xjointlocal = vrotate(xparentjoint + vrotate(constraint.constraints[1].vertices[1], qparentjoint) - xparentbody, inv(qparentbody))
                qjointlocal = qparentbody \ qparentjoint * constraint.constraints[2].qoffset

                # store joint's x and q in world frame
                xjoint = xparentbody + vrotate(xjointlocal, qparentbody)
                qjoint = qparentbody * qjointlocal
                xjointlist[constraint.id] = xjoint
                qjointlist[constraint.id] = qjoint

                # difference to parent body (parentbody)
                qoffset = qjointlocal * qbodylocal

                # actual joint properties
                p1 = xjointlocal # in parent's (parentbody) frame
                p2 = vrotate(-xbodylocal, inv(qbodylocal)) # in body frame (xbodylocal and qbodylocal are both relative to the same (joint) frame -> rotationg by inv(body.q) gives body frame)
                constraint.constraints[1].vertices = (p1, p2)

                V3 = vrotate(constraint.constraints[2].V3', qjointlocal) # in parent's (parentbody) frame
                V12 = (svd(skew(V3)).Vt)[1:2,:]
                constraint.constraints[2].V3 = V3'
                constraint.constraints[2].V12 = V12
                constraint.constraints[2].qoffset = qoffset # in parent's (parentbody) frame

                # actual body properties
                setPosition!(body) # set everything to zero
                setPosition!(parentbody, body, p1 = p1, p2 = p2, Δq = qoffset)
                xbody = body.state.xc
                qbody = body.state.qc


                # shape relative
                if shape !== nothing
                    shape.xoffset = vrotate(xjoint + vrotate(shape.xoffset, qjoint) - xbody, inv(qbody))
                    shape.qoffset = qoffset \ qjointlocal * shape.qoffset
                end
            end
        end

        return mechanism, shapes
    end
end

mutable struct LinearMechanism{T,Nn,Nb,Ne,Ni} <: AbstractMechanism{T,Nn,Nb,Ne,Ni}
    ## Mechanism attributes
    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    graph::Graph{Nn}
    ldu::SparseLDU{T}

    # TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T

    Δt::T
    g::T

    α::T
    μ::T
    
    ## LinearMechanism attributes
    A::AbstractMatrix{T}
    Bu::AbstractMatrix{T}
    Bλ::AbstractMatrix{T}
    G::AbstractMatrix{T}

    xd::Vector{SVector{3,T}}
    vd::Vector{SVector{3,T}}
    qd::Vector{UnitQuaternion{T}}
    ωd::Vector{SVector{3,T}}
    # Fτd::Vector{SVector{3,T}}

    z::Vector{T}
    zsol::Vector{Vector{T}}
    Δz::Vector{T}
    λ::Vector{T}
    λsol::Vector{Vector{T}}
    Δλ::Vector{T}
    u::Vector{T}


    function LinearMechanism(mechanism::Mechanism{T,Nn,Nb,Ne,Ni}, xd, vd, qd, ωd, Fτd, eqcids) where {T,Nn,Nb,Ne,Ni}

        A, Bu, Bλ, G = linearsystem(mechanism, xd, vd, qd, ωd, Fτd, getid.(mechanism.bodies), eqcids)

        z = zeros(T,Nb*12)
        zsol = [zeros(T,Nb*12) for i=1:2]
        Δz = zeros(T,Nb*12)

        nc = 0
        for eqc in mechanism.eqconstraints
            nc += length(eqc)
        end
        λ = zeros(T,nc)
        λsol = [zeros(T,nc) for i=1:2]
        Δλ = zeros(T,nc)
        
        u = zeros(T,size(Bu)[2])

        new{T,Nn,Nb,Ne,Ni}([getfield(mechanism,i) for i=1:getfieldnumber(mechanism)]..., A, Bu, Bλ, G, xd, vd, qd, ωd, z, zsol, Δz, λ, λsol, Δλ, u)
    end

    function LinearMechanism(mechanism::Mechanism{T,Nn,Nb,Ne,Ni}; 
        xd = [mechanism.bodies[i].state.xc for i=1:Nb],
        vd = [mechanism.bodies[i].state.vc for i=1:Nb],
        qd = [mechanism.bodies[i].state.qc for i=1:Nb],
        ωd = [mechanism.bodies[i].state.ωc for i=1:Nb],
        eqcids = [],
        Fτd = []
    ) where {T,Nn,Nb,Ne,Ni}

        return LinearMechanism(mechanism, xd, vd, qd, ωd, Fτd, eqcids)
    end
end


function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, mechanism::AbstractMechanism{T,Nn,Nb,Ne,0}) where {T,Nn,Nb,Ne}
    summary(io, mechanism)
    println(io, " with ", Nb, " bodies and ", Ne, " constraints")
    println(io, " Δt: "*string(mechanism.Δt))
    println(io, " g:  "*string(mechanism.g))
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, mechanism::AbstractMechanism{T,Nn,Nb,Ne,Ni}) where {T,Nn,Nb,Ne,Ni}
    summary(io, mechanism)
    println(io, " with ", Nb, " bodies, ", Ne, " equality constraints, and ", Ni, " inequality constraints")
    println(io, " Δt: "*string(mechanism.Δt))
    println(io, " g:  "*string(mechanism.g))
end
