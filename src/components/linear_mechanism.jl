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

    z0::Vector{T}
    z1::Vector{T}
    Δz::Vector{T}
    λ0::Vector{T}
    λ1::Vector{T}
    Δλ::Vector{T}


    function LinearMechanism(mechanism::Mechanism{T,N,Nb,Ne,Ni}; 
            xd = [szeros(T,3) for i=1:Nb],
            vd = [szeros(T,3) for i=1:Nb],
            qd = [one(UnitQuaternion{T}) for i=1:Nb],
            ωd = [szeros(T,3) for i=1:Nb]
        ) where {T,N,Nb,Ne,Ni}

        A, Bu, Bλ, G = linearsystem(mechanism, xd, vd, qd, ωd, [], getid.(mechanism.bodies), [])

        z0 = zeros(T,Nb*12)
        z1 = zeros(T,Nb*12)
        Δz = zeros(T,Nb*12)

        nc = 0
        for eqc in mechanism.eqconstraints
            nc += length(eqc)
        end
        λ0 = zeros(T,nc)
        λ1 = zeros(T,nc)
        Δλ = zeros(T,nc)


        new{T,N,Nb,Ne,Ni}([getfield(mechanism,i) for i=1:getfieldnumber(mechanism)]..., A, Bu, Bλ, G, xd, vd, qd, ωd, z0, z1, Δz, λ0, λ1, Δλ)
    end
end
