function setPosition!(mechanism::Mechanism{T}, body::Body{T};x::AbstractVector{T} = SVector{3,T}(0, 0, 0),q::Quaternion{T} = Quaternion{T}()) where T
    for i = 1:mechanism.No
        body.state.xd[i] = x
        body.state.qd[i] = q
    end
end

function setPosition!(mechanism::Mechanism{T}, body1::Body{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δx::AbstractVector{T} = SVector{3,T}(0, 0, 0),Δq::Quaternion{T} = Quaternion{T}()) where T

    q2 = body1.state.qd[1] * Δq
    x2 = body1.state.xd[1] + vrotate(SVector{3,T}(p1 + Δx), body1.state.qd[1]) - vrotate(SVector{3,T}(p2), q2)

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
    body.state.xd[1] = body.state.xd[2] - v * Δt
    body.state.qd[1] = Δt/2 * (body.state.qd[2] * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), -SVector{3,T}(ω))) # this accounts for the non-unit-quaternion inverse (q2 * conj(ωbar))
    
    tos!(body, v, ω)
    return
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body1::Body{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δv::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δω::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T

    Δt = mechanism.Δt
    
    v1 = getv1(body1, Δt)
    ω1 = getω1(body1, Δt) # in local coordinates
    
    vp1 = v1 + vrotate(convert(SVector{3,Float64},cross(ω1,p1)),body1.state.qd[2])
    ωp1 = vrotate(ω1,body1.state.qd[2]) # in world coordinates

    vp2 = vp1 + vrotate(convert(SVector{3,Float64},Δv),body1.state.qd[2])
    ωp2 = ωp1 + vrotate(convert(SVector{3,Float64},Δω),body2.state.qd[2]) # in world coordinates

    v2 = vp2 + vrotate(convert(SVector{3,Float64},cross(vrotate(ωp2,inv(body2.state.qd[2])),-p2)),body2.state.qd[2])
    ω2 = vrotate(ωp2,inv(body2.state.qd[2])) # in local coordinates

    setVelocity!(mechanism, body2;v = v2,ω = ω2)
end

# Assumes first order integrator
# TODO higher order integrator
function setVelocity!(mechanism::Mechanism{T}, body1::Origin{T}, body2::Body{T};
    p1::AbstractVector{T} = SVector{3,T}(0, 0, 0), p2::AbstractVector{T} = SVector{3,T}(0, 0, 0), Δv::AbstractVector{T} = SVector{3,T}(0, 0, 0),Δω::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T
    
    vp2 = Δv
    ωp2 = vrotate(convert(SVector{3,Float64},Δω),body2.state.qd[2]) # in world coordinates

    v2 = vp2 + cross(ωp2,-p2)
    ω2 = vrotate(ωp2,inv(body2.state.qd[2])) # in local coordinates

    setVelocity!(mechanism, body2;v = v2,ω = ω2)
end

function setForce!(mechanism::Mechanism{T}, body::Body{T};F::AbstractVector{T} = SVector{3,T}(0, 0, 0),r::AbstractVector{T} = SVector{3,T}(0, 0, 0),τ::AbstractVector{T} = SVector{3,T}(0, 0, 0)) where T
    τ += vrotate(torqueFromForce(F, r),inv(body.state.qd[2])) # in local coordinates
    setForce!(body, F, τ, mechanism.No)
end