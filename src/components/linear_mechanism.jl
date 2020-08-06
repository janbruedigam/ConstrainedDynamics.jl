mutable struct LinearMechanism{T,N,Nb,Ne,Ni} <: AbstractMechanism{T,N,Nb,Ne,Ni}
    ## Mechanism attributes
    origin::Origin{T}
    bodies::UnitDict{Base.OneTo{Int64},Body{T}}
    eqconstraints::UnitDict{UnitRange{Int64},<:EqualityConstraint{T}}
    ineqconstraints::UnitDict{UnitRange{Int64},<:InequalityConstraint{T}}

    graph::Graph{N}
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


    function LinearMechanism(mechanism::Mechanism{T,N,Nb,Ne,Ni}, xd, vd, qd, ωd, Fτd, eqcids) where {T,N,Nb,Ne,Ni}

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

        new{T,N,Nb,Ne,Ni}([getfield(mechanism,i) for i=1:getfieldnumber(mechanism)]..., A, Bu, Bλ, G, xd, vd, qd, ωd, z, zsol, Δz, λ, λsol, Δλ, u)
    end

    function LinearMechanism(mechanism::Mechanism{T,N,Nb,Ne,Ni}; 
        xd = [mechanism.bodies[i].state.xc for i=1:Nb],
        vd = [mechanism.bodies[i].state.vc for i=1:Nb],
        qd = [mechanism.bodies[i].state.qc for i=1:Nb],
        ωd = [mechanism.bodies[i].state.ωc for i=1:Nb],
        eqcids = [],
        Fτd = []
    ) where {T,N,Nb,Ne,Ni}

        return LinearMechanism(mechanism, xd, vd, qd, ωd, Fτd, eqcids)
    end
end
