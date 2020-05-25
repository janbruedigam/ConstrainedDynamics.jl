mutable struct Friction{T} <: Bound{T}
    Nx::Adjoint{T,SVector{6,T}}
    D::SMatrix{2,6,T,12}
    cf::T
    b::SVector{2,T}
    offset::SVector{6,T}


    function Friction(body::Body{T}, normal::AbstractVector{T}, cf::T;offset::AbstractVector{T} = zeros(3)) where T
        @assert cf>0
        normal = normal / norm(normal)

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        A = Array(svd(skew(normal)).V) # gives two plane vectors
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        D = [A[:,1:2];zeros(3,2)]'
        offset = [offset;0;0;0]

        new{T}(Nx, D, cf, zeros(2),offset), body.id
    end
end


@inline function g(friction::Friction, body::Body, Δt, No)
    friction.Nx[SVector(1, 2, 3)]' * (getx2(body, Δt) - friction.offset[SVector(1, 2, 3)])
end

@inline ∂g∂pos(friction::Friction, No) = friction.Nx
@inline ∂g∂vel(friction::Friction, Δt, No) = friction.Nx * Δt

@inline additionalforce(friction::Friction) = friction.D'*friction.b

@inline function schurf(ineqc, friction::Friction, i, body::Body, μ, Δt, No)
    φ = g(friction, body, Δt, No)

    γ1 = ineqc.solγ1[i]
    s1 = ineqc.sols1[i]

    return friction.Nx' * (γ1 / s1 * φ - μ / s1)
end

@inline function schurD(ineqc, friction::Friction, i, body::Body, Δt)
    Nx = friction.Nx
    Nv = Δt * Nx

    γ1 = ineqc.solγ1[i]
    s1 = ineqc.sols1[i]

    return Nx' * γ1 / s1 * Nv
end

# Direct stuff
@inline function calcFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    No = mechanism.No
    cf = friction.cf
    γ1 = ineqc.solγ1[i]
    D = friction.D

    f = body.f
    v = getv2(body)
    ω = getω2(body)
    body.state.vc[2] = @SVector zeros(3)
    body.state.ωc[2] = @SVector zeros(3)
    dyn = dynamics(mechanism, body)
    body.state.vc[2] = v
    body.state.ωc[2] = ω
    body.f = f

    b0 = D*dyn + friction.b # remove old friction force

    if norm(b0) > cf*γ1
        friction.b = b0/norm(b0)*cf*γ1
    else
        friction.b = b0
    end    
    return
end