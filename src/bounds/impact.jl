mutable struct Impact{T} <: Bound{T}
    Nx::Adjoint{T,SVector{6,T}}
    offset::SVector{6,T}


    function Impact(body::Body{T}, normal::AbstractVector{T};offset::AbstractVector{T} = zeros(3)) where T
        normal = normal / norm(normal)

        # Derived from plane equation
        A = Array(svd(skew(normal)).V)
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        offset = [offset;0;0;0]

        new{T}(Nx, offset), body.id
    end
end


@inline function g(impact::Impact, body::Body, dt, No)
    impact.Nx[SVector(1, 2, 3)]' * (getx3(body, dt) - impact.offset[SVector(1, 2, 3)])
end

@inline ∂g∂pos(impact::Impact, No) = impact.Nx
@inline ∂g∂vel(impact::Impact, dt, No) = impact.Nx * dt

@inline function schurf(ineqc, impact::Impact, i, body::Body, μ, dt, No)
    φ = g(impact, body, dt, No)

    γ1 = ineqc.γ1[i]
    s1 = ineqc.s1[i]

    return impact.Nx' * (γ1 / s1 * φ - μ / s1)
end

@inline function schurD(ineqc, impact::Impact, i, body::Body, dt)
    Nx = impact.Nx
    Nv = dt * Nx

    γ1 = ineqc.γ1[i]
    s1 = ineqc.s1[i]

    return Nx' * γ1 / s1 * Nv
end
