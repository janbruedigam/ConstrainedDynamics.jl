# xck = (xdk+1 + xdk)/2
# vck = (xdk+1 - xdk)/Δt

# qck = slerp(qdk+1,qdk,0.5)
# ωck = 2 V inv(qdk) (qdk+1 - qdk)/Δt
# ωckw = sqrt((2/Δt)^2 - ωckᵀωck) - 2/Δt

# Continuous values
@inline getx1(body::Body) = (body.state.xd[2] + body.state.xd[1])/2
@inline getq1(body::Body) = slerp(body.state.qd[2],body.state.qd[1],0.5)
@inline getx2(body::Body, Δt) = body.state.xd[2] + 1/2*getv2(body) * Δt
@inline getq2(body::Body, Δt) = Quaternion(Lmat(body.state.qd[2]) * ωbar(body.state.ωc[2]/2, Δt))
@inline getv1(body::Body) = body.state.vc[1]
@inline getω1(body::Body) = body.state.ωc[1]
@inline getv2(body::Body) = body.state.vc[2]
@inline getω2(body::Body) = body.state.ωc[2]

# Discrete values
@inline getxd3(body::Body, Δt) = body.state.xd[2] + getv2(body) * Δt
@inline getqd3(body::Body, Δt) = Quaternion(Lmat(body.state.qd[2]) * ωbar(body.state.ωc[2], Δt))

@inline function derivωbar(ω::SVector{3,T}, Δt) where T
    msq = -sqrt(4 / Δt^2 - dot(ω, ω))
    Δt / 2 * [ω' / msq; SMatrix{3,3,T,9}(I)]
end

@inline function ωbar(ω, Δt)
    Δt / 2 * Quaternion(sqrt(4 / Δt^2 - dot(ω, ω)), ω)
end

@inline function setForce!(body::Body, F, τ)
    body.F[1] = F
    body.F[2] = F
    body.τ[1] = τ
    body.τ[2] = τ
end


@inline function discretizestate!(body::Body, Δt)
    state = body.state
    xc = state.xc[1]
    qc = state.qc[1]
    vc = state.vc[1]
    ωc = state.ωc[1]

    state.xd[1] = xc - 1/2*vc*Δt
    state.xd[2] = xc + 1/2*vc*Δt
    state.qd[1] = Quaternion(Lmat(qc) * ωbar(-ωc/2, Δt))
    state.qd[2] = Quaternion(Lmat(qc) * ωbar(ωc/2, Δt))

    # just to set to some reasonable value
    state.xc[2] = xc + vc*Δt
    state.qc[2] = Quaternion(Lmat(qc) * ωbar(ωc, Δt))
    state.vc[2] = vc
    state.ωc[2] = ωc

    return
end

@inline function dynamics(mechanism, body::Body{T}) where T
    Δt = mechanism.Δt

    ezg = SVector{3,T}(0, 0, -mechanism.g)
    dynT = body.m * ((getv2(body) - getv1(body)) / Δt + ezg) - 1/2*(body.F[1]+body.F[2])

    J = body.J
    ω1 = getω1(body)
    ω2 = getω2(body)
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    dynR = skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1) - (body.τ[1] + body.τ[2])

    body.f = [dynT;dynR]

    for cid in connections(mechanism.graph, body.id)
        GtλTof!(mechanism, body, geteqconstraint(mechanism, cid))
    end

    for cid in ineqchildren(mechanism.graph, body.id)
        NtγTof!(mechanism, body, getineqconstraint(mechanism, cid))
    end

    return body.f
end

@inline function ∂dyn∂pos(body::Body{T}, Δt) where T
    J = body.J
    ω2 = getω2(body)
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    dynT = SMatrix{3,3,T,9}(2 * body.m / Δt^2 * I)
    dynR = 2 * (skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2)) * 2/Δt * VLᵀmat(body.state.qd[2])*LVᵀmat(getq2(body,Δt))

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end

@inline function ∂dyn∂vel(body::Body{T}, Δt) where T
    J = body.J
    ω2 = getω2(body)
    sq = sqrt(4 / Δt^2 - ω2' * ω2)

    dynT = SMatrix{3,3,T,9}(body.m / Δt * I)
    dynR = skewplusdiag(ω2, sq) * J - J * ω2 * (ω2' / sq) - skew(J * ω2)

    Z = @SMatrix zeros(T, 3, 3)

    return [[dynT; Z] [Z; dynR]]
end



####
# Tests
####

@inline function discretizestate!(body,x1,q1,v1,v2,ω1,ω2,Δt)
    state = body.state

    state.xc[1] = x1
    state.qc[1] = q1
    state.vc[1] = v1
    state.ωc[1] = ω1
    state.xc[2] = x1 + 1/2*(v1+v2)*Δt
    state.qc[2] = q1 * ωbar(ω1/2,Δt) * ωbar(ω2/2,Δt)
    state.vc[2] = v2
    state.ωc[2] = ω2

    state.xd[1] = x1 - 1/2*v1*Δt
    state.xd[2] = x1 + 1/2*v1*Δt
    state.qd[1] = q1 / ωbar(ω1/2,Δt)
    state.qd[2] = q1 * ωbar(ω1/2,Δt)

    return
end

function velTpos(vars)
    ezg = SVector{3,Float64}(0, 0, 9.81)
    Δt = 0.01
    m = 1.

    xc1 = vars[1:3]
    vc2 = vars[4:6]

    xd2 = xc1
    xd3 = xd2 + vc2 * Δt
    
    2/Δt
end

function dynTvel(vars)
    ezg = SVector{3,Float64}(0, 0, 9.81)
    Δt = 0.01
    m = 1.

    vc2 = vars[1:3]
    vc1 = vars[4:6]

    m * ((vc2 - vc1) / Δt + ezg)
end

function velRpos(vars)
    Δt = 0.01

    qc1 = Quaternion(SVector(vars[1:4]...))
    ωc2 = vars[5:7]

    qd2 = qc1 * ωbar(ωc2,Δt)
    qc2 = qd2 * ωbar(ωc2/2,Δt)

    2 * 2/Δt * VLᵀmat(qd2) * Lmat(Quaternion(qc2)) * Vᵀmat()
end

function dynRvel(vars)
    Δt = 0.01

    J = diagm(ones(3))
    ωc2 = vars[1:3]
    ωc1 = vars[4:6]
    sq2 = sqrt(4 / Δt^2 - ωc2' * ωc2)
    sq1 = sqrt(4 / Δt^2 - ωc1' * ωc1)
    
    skewplusdiag(ωc2, sq2) * (J * ωc2) - skewplusdiag(ωc1, sq1) * (J * ωc1)
end
