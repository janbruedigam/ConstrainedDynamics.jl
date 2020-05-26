@inline function discretizestate!(body,x1,q1,v1,v2,ω1,ω2,Δt)
    state = body.state

    state.xc[1] = x1
    state.qc[1] = q1
    state.vc[1] = v1
    state.ωc[1] = ω1
    state.vc[2] = v2
    state.ωc[2] = ω2

    state.xd[1] = x1 - 1/2*v1*Δt
    state.xd[2] = x1 + 1/2*v1*Δt
    state.qd[1] = q1 / ωbar(ω1/2,Δt)
    state.qd[2] = q1 * ωbar(ω1/2,Δt)

    return
end



# Dyn diff test

function dynTvel(vars)
    ezg = SVector{3,Float64}(0, 0, 9.81)
    Δt = 0.01
    m = 1.

    vc1 = vars[1:3]
    vc2 = vars[4:6]

    m * ((vc2 - vc1) / Δt + ezg)
end

function dynRvel(vars)
    Δt = 0.01

    J = diagm(ones(3))
    ωc1 = vars[1:3]
    ωc2 = vars[4:6]
    sq2 = sqrt(4 / Δt^2 - ωc2' * ωc2)
    sq1 = sqrt(4 / Δt^2 - ωc1' * ωc1)
    
    skewplusdiag(ωc2, sq2) * (J * ωc2) - skewplusdiag(ωc1, sq1) * (J * ωc1)
end


# Joint diff test

function transfunc3pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc2pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    V1 = vars[21:23]
    V2 = vars[24:26]
    V12 = [V1';V2']

    V12 * vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc1pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    V3 = vars[21:23]

    V3' * vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end


function transfunc3vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = Quaternion(SVector(vars[4:7]...))
    va2 = vars[8:10]
    ωa2 = vars[11:13]

    xb1 = vars[14:16]
    qb1 = Quaternion(SVector(vars[17:20]...))
    vb2 = vars[21:23]
    ωb2 = vars[24:26]

    xa2 = xa1 + 1/2 * va2 * Δt
    qa2 = Quaternion(Lmat(qa1) * ωbar(1/2 * ωa2, Δt))

    xb2 = xb1 + 1/2 * vb2 * Δt
    qb2 = Quaternion(Lmat(qb1) * ωbar(1/2 * ωb2, Δt))

    pa = vars[27:29]
    pb = vars[30:32]

    vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end

function transfunc2vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = Quaternion(SVector(vars[4:7]...))
    va2 = vars[8:10]
    ωa2 = vars[11:13]

    xb1 = vars[14:16]
    qb1 = Quaternion(SVector(vars[17:20]...))
    vb2 = vars[21:23]
    ωb2 = vars[24:26]

    xa2 = xa1 + 1/2 * va2 * Δt
    qa2 = Quaternion(Lmat(qa1) * ωbar(1/2 * ωa2, Δt))

    xb2 = xb1 + 1/2 * vb2 * Δt
    qb2 = Quaternion(Lmat(qb1) * ωbar(1/2 * ωb2, Δt))

    pa = vars[27:29]
    pb = vars[30:32]

    V1 = vars[33:35]
    V2 = vars[36:38]
    V12 = [V1';V2']

    V12 * vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end

function transfunc1vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = Quaternion(SVector(vars[4:7]...))
    va2 = vars[8:10]
    ωa2 = vars[11:13]

    xb1 = vars[14:16]
    qb1 = Quaternion(SVector(vars[17:20]...))
    vb2 = vars[21:23]
    ωb2 = vars[24:26]

    xa2 = xa1 + 1/2 * va2 * Δt
    qa2 = Quaternion(Lmat(qa1) * ωbar(1/2 * ωa2, Δt))

    xb2 = xb1 + 1/2 * vb2 * Δt
    qb2 = Quaternion(Lmat(qb1) * ωbar(1/2 * ωb2, Δt))

    pa = vars[27:29]
    pb = vars[30:32]

    V3 = vars[33:35]

    V3' * vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end



function rotfunc3pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    offset = Quaternion(SVector(vars[15:18]...))

    VLᵀmat(offset) * Lᵀmat(qa) * qb
end

function rotfunc2pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    offset = Quaternion(SVector(vars[15:18]...))

    V1 = vars[19:21]
    V2 = vars[22:24]
    V12 = [V1';V2']

    V12 * VLᵀmat(offset) * Lᵀmat(qa) * qb
end

function rotfunc1pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    offset = Quaternion(SVector(vars[15:18]...))

    V3 = vars[19:21]

    V3' * VLᵀmat(offset) * Lᵀmat(qa) * qb
end


function rotfunc3vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = Quaternion(SVector(vars[4:7]...))
    va2 = vars[8:10]
    ωa2 = vars[11:13]

    xb1 = vars[14:16]
    qb1 = Quaternion(SVector(vars[17:20]...))
    vb2 = vars[21:23]
    ωb2 = vars[24:26]

    xa2 = xa1 + 1/2 * va2 * Δt
    qa2 = Quaternion(Lmat(qa1) * ωbar(1/2 * ωa2, Δt))

    xb2 = xb1 + 1/2 * vb2 * Δt
    qb2 = Quaternion(Lmat(qb1) * ωbar(1/2 * ωb2, Δt))

    offset = Quaternion(SVector(vars[27:30]...))

    VLᵀmat(offset) * Lᵀmat(qa2) * qb2
end

function rotfunc2vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = Quaternion(SVector(vars[4:7]...))
    va2 = vars[8:10]
    ωa2 = vars[11:13]

    xb1 = vars[14:16]
    qb1 = Quaternion(SVector(vars[17:20]...))
    vb2 = vars[21:23]
    ωb2 = vars[24:26]

    xa2 = xa1 + 1/2 * va2 * Δt
    qa2 = Quaternion(Lmat(qa1) * ωbar(1/2 * ωa2, Δt))

    xb2 = xb1 + 1/2 * vb2 * Δt
    qb2 = Quaternion(Lmat(qb1) * ωbar(1/2 * ωb2, Δt))

    offset = Quaternion(SVector(vars[27:30]...))

    V1 = vars[31:33]
    V2 = vars[34:36]
    V12 = [V1';V2']    

    V12 * VLᵀmat(offset) * Lᵀmat(qa2) * qb2
end

function rotfunc1vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = Quaternion(SVector(vars[4:7]...))
    va2 = vars[8:10]
    ωa2 = vars[11:13]

    xb1 = vars[14:16]
    qb1 = Quaternion(SVector(vars[17:20]...))
    vb2 = vars[21:23]
    ωb2 = vars[24:26]

    xa2 = xa1 + 1/2 * va2 * Δt
    qa2 = Quaternion(Lmat(qa1) * ωbar(1/2 * ωa2, Δt))

    xb2 = xb1 + 1/2 * vb2 * Δt
    qb2 = Quaternion(Lmat(qb1) * ωbar(1/2 * ωb2, Δt))

    offset = Quaternion(SVector(vars[27:30]...))

    V3 = vars[31:33]

    V3' * VLᵀmat(offset) * Lᵀmat(qa2) * qb2
end