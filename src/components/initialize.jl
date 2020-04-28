function setPosition!(mechanism::Mechanism{T}, body::Body{T};x::AbstractVector{T} = SVector{3,T}(0, 0, 0),q::Quaternion{T} = Quaternion{T}()) where T
    for i = 1:mechanism.No
        body.x[i] = x
        body.q[i] = q
    end
end

function setPosition!(mechanism::Mechanism{T}, body1::Body{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δx::AbstractVector{T} = SVector{3,T}(0, 0, 0),Δq::Quaternion{T} = Quaternion{T}()) where T

    q2 = body1.q[1] * Δq
    x2 = body1.x[1] + vrotate(SVector{3,T}(p1 + Δx), body1.q[1]) - vrotate(SVector{3,T}(p2), q2)

    setPosition!(mechanism, body2;x = x2,q = q2)
end

function setPosition!(mechanism::Mechanism{T}, body1::Origin{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δx::AbstractVector{T} = SVector{3,T}(0, 0, 0),Δq::Quaternion{T} = Quaternion{T}()) where T

    q2 = Δq
    x2 = p1 + Δx - vrotate(SVector{3,T}(p2), q2)


    setPosition!(mechanism, body2;x = x2,q = q2)
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body::Body{T};v::AbstractVector{T} = SVector{3,T}(0, 0, 0),ω::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T
    Δt = mechanism.Δt
    body.x[1] = body.x[2] - v * Δt
    body.q[1] = Δt/2 * (body.q[2] * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), -SVector{3,T}(ω))) # this accounts for the non-unit-quaternion inverse (q2 * conj(ωbar))
    body.s0 = [v;ω]
    s0tos1!(body)
    return
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body1::Body{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δv::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δω::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T

    Δt = mechanism.Δt
    
    v1 = getv1(body1, Δt)
    ω1 = getω1(body1, Δt)
    
    vp1 = v1 + cross(ω1,p1)
    ωp1 = ω1

    vp2 = vp1 + Δv
    ωp2 = ωp1 + Δω

    v2 = vp2 + cross(ωp2,-p2)
    ω2 = ωp2

    setVelocity!(mechanism, body2;v = v2,ω = ω2)
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body1::Origin{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δv::AbstractVector{T} = SVector{3,T}(0, 0, 0),Δω::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T
    
    vp2 = Δv
    ωp2 = Δω

    v2 = vp2 + cross(ωp2,-p2)
    ω2 = ωp2

    setVelocity!(mechanism, body2;v = v2,ω = ω2)
end

function setForce!(mechanism::Mechanism{T}, body::Body{T};F::AbstractVector{T} = SVector{3,T}(0, 0, 0),r::AbstractVector{T} = SVector{3,T}(0, 0, 0),τ::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T
    τ += torqueFromForce(F, r)
    setForce!(body, F, τ, mechanism.No)
end