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
    friction.Nx[SVector(1, 2, 3)]' * (getx3(body, Δt) - friction.offset[SVector(1, 2, 3)])
end

@inline ∂g∂pos(friction::Friction, No) = friction.Nx
@inline ∂g∂vel(friction::Friction, Δt, No) = friction.Nx * Δt

@inline function schurf(ineqc, friction::Friction, i, body::Body, μ, Δt, No)
    φ = g(friction, body, Δt, No)

    γ1 = ineqc.γ1[i]
    s1 = ineqc.s1[i]

    return friction.Nx' * (γ1 / s1 * φ - μ / s1)
end

@inline function schurD(ineqc, friction::Friction, i, body::Body, Δt)
    Nx = friction.Nx
    Nv = Δt * Nx

    γ1 = ineqc.γ1[i]
    s1 = ineqc.s1[i]

    return Nx' * γ1 / s1 * Nv
end

# Smooth stuff
@inline function setFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
    Δt = mechanism.Δt
    No = mechanism.No
    M = getM(body)
    v = body.s1
    cf = friction.cf
    γ1 = ineqc.γ1[i]
    D = friction.D

    B = D'*friction.b
    F = body.F[No] - B[SVector(1,2,3)]
    τ = body.τ[No] - B[SVector(4,5,6)]
    setForce!(body,F,τ,No)

    ψ = Δt*norm(D*v)
    
    f = body.f
    body.s1 = @SVector zeros(6)
    dyn = D/M*dynamics(body,mechanism)*Δt^2
    body.s1 = v
    body.f = f
    
    X = D/M*D' * Δt^2 + I*(ψ/(cf*γ1))

    friction.b = X\dyn
    B = D'*friction.b
    F += B[SVector(1,2,3)]
    τ += B[SVector(4,5,6)]
    setForce!(body,F,τ,No)
    return
end

# Prox stuff
# @inline function setFrictionForce!(mechanism, ineqc, friction::Friction, i, body::Body)
#     Δt = mechanism.Δt
#     No = mechanism.No
#     M = getM(body)
#     v = body.s1
#     cf = friction.cf
#     γ1 = ineqc.γ1[i]
#     D = friction.D

#     B = D'*friction.b
#     F = body.F[No] - B[SVector(1,2,3)]
#     τ = body.τ[No] - B[SVector(4,5,6)]
#     setForce!(body,F,τ,No)

#     ψ = Δt*norm(D*v)
    
#     f = body.f
#     body.s1 = @SVector zeros(6)
#     dyn = dynamics(body,mechanism)
#     # body.s1 = v
#     # body.f = f

#     friction.b = D*dyn
#     if norm(friction.b) > cf*γ1
#         friction.b = friction.b/norm(friction.b)*cf*γ1
#     end

#     B = D'*friction.b
#     F += B[SVector(1,2,3)]
#     τ += B[SVector(4,5,6)]
#     setForce!(body,F,τ,No)

#     body.s1 = v
#     body.f = f
#     # v = M\dynamics(body,mechanism)*Δt
#     # body.s1 = v
#     # dynamics(body,mechanism)

#     return
# end