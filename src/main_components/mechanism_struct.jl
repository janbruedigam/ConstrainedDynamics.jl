abstract type AbstractMechanism{T,Nn,Nb,Ne,Ni} end

"""
$(TYPEDEF)

A `Mechanism` contains the [`Origin`](@ref), [`Body`](@ref)s, and [`EqualityConstraint`](@ref)s of a system and can be used for simulation.
# Important attributes
* `origin`:        The origin of a mechanism.
* `bodies`:        The bodies of a mechanism (Dict).
* `eqconstraints`: The equality constraints (joints) of a mechanism (Dict).
* `Δt`:            The time step of the mechanism.
* `g`:             The gravitational acceleration in z-direction.

# Constuctors
    Mechanism(origin, bodies; Δt, g)
    Mechanism(origin, bodies, eqcs; Δt, g)
    Mechanism(urdf_filename; floating, Δt, g)
"""
mutable struct Mechanism{T,Nn,Nb,Ne,Ni} <: AbstractMechanism{T,Nn,Nb,Ne,Ni}
    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    graph::Graph{Nn}
    structure::Structure

    # TODO remove once EqualityConstraint is homogenous
    normf::T
    normΔs::T

    Δt::T
    g::T

    α::T
    μ::T


    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},
            eqcs::Vector{<:EqualityConstraint{T}}, ineqcs::Vector{<:InequalityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81
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

        if Nb < Ne
            @info "More constraints than bodies. Potentially bad behavior."
        end


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

        structure = Structure(origin, bodies, eqcs, ineqcs)
        graph = Graph(origin, bodies, eqcs, ineqcs)

        bodies = UnitDict(bodies)
        eqcs = UnitDict((eqcs[1].id):(eqcs[Ne].id), eqcs)
        if Ni > 0
            ineqcs = UnitDict((ineqcs[1].id):(ineqcs[Ni].id), ineqcs)
        else
            ineqcs = UnitDict(0:0, ineqcs)
        end

        α = 1
        μ = 1

        new{T,Nn,Nb,Ne,Ni}(origin, bodies, eqcs, ineqcs, graph, structure, normf, normΔs, Δt, g, α, μ)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},eqcs::Vector{<:EqualityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        ineqcs = InequalityConstraint{T}[]
        return Mechanism(origin, bodies, eqcs, ineqcs, Δt = Δt, g = g)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}},ineqcs::Vector{<:InequalityConstraint{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(Floating(origin, body)))
        end
        return Mechanism(origin, bodies, eqc, ineqcs, Δt = Δt, g = g)
    end

    function Mechanism(origin::Origin{T},bodies::Vector{Body{T}};
            Δt::Real = .01, g::Real = -9.81
        ) where T

        eqc = EqualityConstraint{T}[]
        for body in bodies
            push!(eqc, EqualityConstraint(Floating(origin, body)))
        end
        return Mechanism(origin, bodies, eqc, Δt = Δt, g = g)
    end

    function Mechanism(filename::AbstractString; floating::Bool=false, type::Type{T} = Float64, Δt::Real = .01, g::Real = -9.81) where T
        origin, links, joints, loopjoints = parse_urdf(filename, floating, T)
        mechanism = Mechanism(origin, links, [joints;loopjoints], Δt = Δt, g = g)
        set_parsed_values!(mechanism, loopjoints)

        return mechanism
    end
end

mutable struct LinearMechanism{T,Nn,Nb,Ne,Ni} <: AbstractMechanism{T,Nn,Nb,Ne,Ni}
    ## Mechanism attributes
    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    graph::Graph{Nn}
    structure::Structure

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
