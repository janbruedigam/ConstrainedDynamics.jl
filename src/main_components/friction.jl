mutable struct Friction{T,N} <: Component{T}
    id::Int64
    name::String

    Bx::SMatrix{4,3,T,12}
    p::SVector{3,T}

    parentid::Union{Int64,Nothing}
    childids::SVector{2,Int64}

    βsol::Vector{SVector{N,T}}

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

        new{T,N}(frictionid, name, Bx, p, body.id, [frictionbound.id; betabound.id], βsol), [impact; frictionbound; betabound]
    end
end


Base.length(::Friction{T,N}) where {T,N} = N

@inline ∂g∂ʳself(mechanism, fric::Friction{T,N}) where {T,N} = szeros(T,N,N)

@inline Bq(Bxmat, p, q) = Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p)
@inline Bmat(Bxmat, p, q) = [Bxmat Bq(Bxmat, p, q)]

function g(mechanism, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    frictionbound = getineqconstraint(mechanism, fric.childids[1])
    betabound = getineqconstraint(mechanism, fric.childids[2])

    x, v, q, ω = fullargssol(body.state)
    ψ = frictionbound.γsol[2]
    η = betabound.γsol[2]

    Bxmat = fric.Bx
    Bqmat = Bq(Bxmat, fric.p, q)
    return Bxmat*v + Bqmat*ω - η .+ ψ
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
