@inline function discretizestate!(body,x1,q1,v1,v2,ω1,ω2,Δt)
    state = body.state

    state.xc[1] = x1
    state.qc[1] = q1
    state.vc[1] = v1
    state.ωc[1] = ω1
    state.xc[2] = x1 + v1*Δt
    state.qc[2] = q1 * ωbar(ω1,Δt)
    state.vc[2] = v2
    state.ωc[2] = ω2

    state.xd[1] = x1 - v1*Δt
    state.xd[2] = x1
    state.qd[1] = q1 / ωbar(ω1,Δt)
    state.qd[2] = q1

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
    
    1/Δt
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

    qd2 = qc1
    qc2 = qd2 * ωbar(ωc2,Δt)

    2/Δt * VLᵀmat(qd2) * Lmat(Quaternion(qc2)) * Vᵀmat()
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