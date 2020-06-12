@inline function discretizestate!(body,x1,q1,v1,v2,ω1,ω2,Δt)
    state = body.state

    state.xc = x1
    state.qc = q1
    state.vc = v1
    state.ωc = ω1
    state.vsol[2] = v2
    state.ωsol[2] = ω2

    state.xk[1] = x1
    state.xk[2] = x1 + v1*Δt
    state.qk[1] = q1
    state.qk[2] = q1 * ωbar(ω1,Δt)

    return
end

@inline getxqkvector(state) = [state.xk[2];params(state.qk[2])]
@inline getxk(state) = state.xk[2]
@inline getqk(state) = state.qk[2]
@inline getstateandvestimate(state) = [getxqkvector(state);state.vc;state.ωc;state.vsol[2];state.ωsol[2]]


# Dyn diff test

function dynTvel(vars)
    ezg = SVector{3,Float64}(0, 0, 9.81)
    Δt = 0.01
    m = 1.

    vc1 = vars[1:3]
    vc2 = vars[4:6]

    return m * ((vc2 - vc1) / Δt + ezg)
end

function dynRvel(vars)
    Δt = 0.01

    J = diagm(ones(3))
    ωc1 = vars[1:3]
    ωc2 = vars[4:6]
    sq2 = sqrt(4 / Δt^2 - ωc2' * ωc2)
    sq1 = sqrt(4 / Δt^2 - ωc1' * ωc1)
    
    return skewplusdiag(ωc2, sq2) * (J * ωc2) - skewplusdiag(ωc1, sq1) * (J * ωc1)
end


# Joint diff test

function transfunc3pos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = UnitQuaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    return vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc2pos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = UnitQuaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    V1 = vars[21:23]
    V2 = vars[24:26]
    V12 = [V1';V2']

    return V12 * vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc1pos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = UnitQuaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    V3 = vars[21:23]

    return V3' * vrotate(xb + vrotate(pb, qb) - (xa + vrotate(pa, qa)), inv(qa))
end


function transfunc3vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = UnitQuaternion(SVector(vars[4:7]...))
    va1 = vars[8:10]
    ωa1 = vars[11:13]
    va2 = vars[14:16]
    ωa2 = vars[17:19]

    xb1 = vars[20:22]
    qb1 = UnitQuaternion(SVector(vars[23:26]...))
    vb1 = vars[27:29]
    ωb1 = vars[30:32]
    vb2 = vars[33:35]
    ωb2 = vars[36:38]

    xa2 = xa1 + va2 * Δt
    qa2 = qa1 * ωbar(ωa2, Δt)

    xb2 = xb1 + vb2 * Δt
    qb2 = qb1 * ωbar(ωb2, Δt)

    pa = vars[39:41]
    pb = vars[42:44]

    return vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end

function transfunc2vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = UnitQuaternion(SVector(vars[4:7]...))
    va1 = vars[8:10]
    ωa1 = vars[11:13]
    va2 = vars[14:16]
    ωa2 = vars[17:19]

    xb1 = vars[20:22]
    qb1 = UnitQuaternion(SVector(vars[23:26]...))
    vb1 = vars[27:29]
    ωb1 = vars[30:32]
    vb2 = vars[33:35]
    ωb2 = vars[36:38]

    xa2 = xa1 + va2 * Δt
    qa2 = qa1 * ωbar(ωa2, Δt)

    xb2 = xb1 + vb2 * Δt
    qb2 = qb1 * ωbar(ωb2, Δt)

    pa = vars[39:41]
    pb = vars[42:44]

    V1 = vars[45:47]
    V2 = vars[48:50]
    V12 = [V1';V2']

    return V12 * vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end

function transfunc1vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = UnitQuaternion(SVector(vars[4:7]...))
    va1 = vars[8:10]
    ωa1 = vars[11:13]
    va2 = vars[14:16]
    ωa2 = vars[17:19]

    xb1 = vars[20:22]
    qb1 = UnitQuaternion(SVector(vars[23:26]...))
    vb1 = vars[27:29]
    ωb1 = vars[30:32]
    vb2 = vars[33:35]
    ωb2 = vars[36:38]

    xa2 = xa1 + va2 * Δt
    qa2 = qa1 * ωbar(ωa2, Δt)

    xb2 = xb1 + vb2 * Δt
    qb2 = qb1 * ωbar(ωb2, Δt)

    pa = vars[39:41]
    pb = vars[42:44]

    V3 = vars[45:47]

    return V3' * vrotate(xb2 + vrotate(pb, qb2) - (xa2 + vrotate(pa, qa2)), inv(qa2))
end



function rotfunc3pos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = UnitQuaternion(SVector(vars[11:14]...))

    offset = UnitQuaternion(SVector(vars[15:18]...))

    return VLᵀmat(offset) * Lᵀmat(qa) * params(qb)
end

function rotfunc2pos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = UnitQuaternion(SVector(vars[11:14]...))

    offset = UnitQuaternion(SVector(vars[15:18]...))

    V1 = vars[19:21]
    V2 = vars[22:24]
    V12 = [V1';V2']

    return V12 * VLᵀmat(offset) * Lᵀmat(qa) * params(qb)
end

function rotfunc1pos(vars)
    xa = vars[1:3]
    qa = UnitQuaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = UnitQuaternion(SVector(vars[11:14]...))

    offset = UnitQuaternion(SVector(vars[15:18]...))

    V3 = vars[19:21]

    return V3' * VLᵀmat(offset) * Lᵀmat(qa) * params(qb)
end


function rotfunc3vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = UnitQuaternion(SVector(vars[4:7]...))
    va1 = vars[8:10]
    ωa1 = vars[11:13]
    va2 = vars[14:16]
    ωa2 = vars[17:19]

    xb1 = vars[20:22]
    qb1 = UnitQuaternion(SVector(vars[23:26]...))
    vb1 = vars[27:29]
    ωb1 = vars[30:32]
    vb2 = vars[33:35]
    ωb2 = vars[36:38]

    xa2 = xa1 + va2 * Δt
    qa2 = qa1 * ωbar(ωa2, Δt)

    xb2 = xb1 + vb2 * Δt
    qb2 = qb1 * ωbar(ωb2, Δt)

    offset = UnitQuaternion(SVector(vars[39:42]...))

    return VLᵀmat(offset) * Lᵀmat(qa2) * params(qb2)
end

function rotfunc2vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = UnitQuaternion(SVector(vars[4:7]...))
    va1 = vars[8:10]
    ωa1 = vars[11:13]
    va2 = vars[14:16]
    ωa2 = vars[17:19]

    xb1 = vars[20:22]
    qb1 = UnitQuaternion(SVector(vars[23:26]...))
    vb1 = vars[27:29]
    ωb1 = vars[30:32]
    vb2 = vars[33:35]
    ωb2 = vars[36:38]

    xa2 = xa1 + va2 * Δt
    qa2 = qa1 * ωbar(ωa2, Δt)

    xb2 = xb1 + vb2 * Δt
    qb2 = qb1 * ωbar(ωb2, Δt)

    offset = UnitQuaternion(SVector(vars[39:42]...))

    V1 = vars[43:45]
    V2 = vars[46:48]
    V12 = [V1';V2']    

    return V12 * VLᵀmat(offset) * Lᵀmat(qa2) * params(qb2)
end

function rotfunc1vel(vars)
    Δt = 0.01

    xa1 = vars[1:3]
    qa1 = UnitQuaternion(SVector(vars[4:7]...))
    va1 = vars[8:10]
    ωa1 = vars[11:13]
    va2 = vars[14:16]
    ωa2 = vars[17:19]

    xb1 = vars[20:22]
    qb1 = UnitQuaternion(SVector(vars[23:26]...))
    vb1 = vars[27:29]
    ωb1 = vars[30:32]
    vb2 = vars[33:35]
    ωb2 = vars[36:38]

    xa2 = xa1 + va2 * Δt
    qa2 = qa1 * ωbar(ωa2, Δt)

    xb2 = xb1 + vb2 * Δt
    qb2 = qb1 * ωbar(ωb2, Δt)

    offset = UnitQuaternion(SVector(vars[39:42]...))

    V3 = vars[43:45]

    return V3' * VLᵀmat(offset) * Lᵀmat(qa2) * params(qb2)
end
