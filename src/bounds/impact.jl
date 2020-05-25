mutable struct Impact{T} <: Bound{T}
    Nx::Adjoint{T,SVector{6,T}}
    offset::SVector{6,T}


    function Impact(body::Body{T}, normal::AbstractVector{T};offset::AbstractVector{T} = zeros(3)) where T
        normal = normal / norm(normal)

        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        A = Array(svd(skew(normal)).V) # gives two plane vectors
        A[:,3] = normal # to ensure correct sign
        Ainv = inv(A)
        ainv3 = Ainv[3,:]
        Nx = [ainv3;0;0;0]'
        offset = [offset;0;0;0]

        new{T}(Nx, offset), body.id
    end
end


@inline function g(impact::Impact, body::Body, Δt, No)
    impact.Nx[SVector(1, 2, 3)]' * (getx2(body, Δt) - impact.offset[SVector(1, 2, 3)])
end

@inline ∂g∂pos(impact::Impact, No) = impact.Nx
@inline ∂g∂vel(impact::Impact, Δt, No) = impact.Nx * Δt

@inline function schurf(ineqc, impact::Impact, i, body::Body, μ, Δt, No)
    φ = g(impact, body, Δt, No)

    γ1 = ineqc.solγ1[i]
    s1 = ineqc.sols1[i]

    return impact.Nx' * (γ1 / s1 * φ - μ / s1)
end

@inline function schurD(ineqc, impact::Impact, i, body::Body, Δt)
    Nx = impact.Nx
    Nv = Δt * Nx

    γ1 = ineqc.solγ1[i]
    s1 = ineqc.sols1[i]

    return Nx' * γ1 / s1 * Nv
end
