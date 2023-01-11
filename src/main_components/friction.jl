mutable struct Friction{T} <: Component{T}
    id::Int64
    name::String

    Bx::SMatrix{4,3,T,12}
    p::SVector{3,T}

    parentid::Int64
    childids::SVector{2,Int64}

    βsol::Vector{SVector{4,T}}
    γsolref::Vector{SVector{1,T}} #TODO this is a reference to the associated impact γ to avoid allocations

    d::SVector{4,T}

    function Friction(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3), name::String="") where T
        N = 4

        Bx = SA{T}[
            1 0 0
            -1 0 0
            0 1 0
            0 -1 0
        ]

        frictionid = getGlobalID()
        impact = InequalityConstraint(Impact(body, normal; p = p, offset = offset))
        frictionbound = InequalityConstraint(FrictionBound(cf, frictionid, impact.id))
        betabound = InequalityConstraint(BetaBound(frictionid))

        βsol = [szeros(T, N) for i=1:2]
        γsolref = impact.γsol

        d = szeros(T, 4)

        new{T}(frictionid, name, Bx, p, body.id, [frictionbound.id; betabound.id], βsol, γsolref, d), [impact; frictionbound; betabound]
    end
end

function resetVars!(fric::Friction{T}) where {T}
    fric.βsol[1] = szeros(T, 4)
    fric.βsol[2] = szeros(T, 4)

    return 
end

Base.length(::Friction) = 4

@inline ∂g∂ʳself(mechanism, fric::Friction{T}) where {T} = szeros(T,4,4)

@inline Bq(Bxmat, p, q) = Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p)
@inline Bmat(Bxmat, p, q) = [Bxmat Bq(Bxmat, p, q)]

function g(mechanism, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    x, v, q, ω = fullargssol(body.state)
    Bxmat = fric.Bx
    Bqmat = Bq(Bxmat, fric.p, q)

    fric.d = Bxmat*v + Bqmat*ω

    for childid in fric.childids
        constraintForceMapping!(mechanism, fric, getineqconstraint(mechanism, childid))
    end
    
    return fric.d
end

@inline function constraintForceMapping!(mechanism, body::Body, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    x, q = posargsk(body.state)

    body.state.d -= Bmat(fric.Bx, fric.p, q)' * fric.βsol[2]
    return
end

@inline function ∂gab∂ʳba(mechanism, body::Body, fric::Friction)
    x, q = posargsk(body.state)
    B = Bmat(fric.Bx, fric.p, q)

    return -B', B
end

@inline function ∂gab∂ʳba(mechanism, fric::Friction, ineqc::InequalityConstraint)
    G = ∂g∂beta(fric, ineqc.constraints[1])

    return G, -G'
end
