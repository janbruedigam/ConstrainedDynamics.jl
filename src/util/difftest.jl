using ForwardDiff
using Rotations
using StaticArrays
using LinearAlgebra

using ConstrainedDynamics
using ConstrainedDynamics: vrotate, Lmat, Vᵀmat, Lᵀmat, VLᵀmat, ∂g∂pos, ∂g∂vel, skew

function transfunc3pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    (xb + vrotate(SVector{3}(pb...), qb)) - (xa + vrotate(SVector{3}(pa), qa))
end

function transfunc3vel(vars)
    Δt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    pa = vars[27:29]
    pb = vars[30:32]

    xa = xa + Δt * va
    xb = xb + Δt * vb
    qa = Quaternion(Δt / 2 * (qa * Quaternion(SVector(sqrt((2 / Δt)^2 - wa' * wa), wa...))))
    qb = Quaternion(Δt / 2 * (qb * Quaternion(SVector(sqrt((2 / Δt)^2 - wb' * wb), wb...))))

    (xb + vrotate(SVector{3}(pb...), qb)) - (xa + vrotate(SVector{3}(pa...), qa))
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

    V12 * vrotate(SVector{3}(((xb + vrotate(SVector{3}(pb...), qb)) - (xa + vrotate(SVector{3}(pa...), qa)))...), inv(qa))
end

function transfunc2vel(vars)
    Δt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    pa = vars[27:29]
    pb = vars[30:32]

    V1 = vars[33:35]
    V2 = vars[36:38]
    V12 = [V1';V2']

    xa = xa + Δt * va
    xb = xb + Δt * vb
    qa = Quaternion(Δt / 2 * (qa * Quaternion(SVector(sqrt((2 / Δt)^2 - wa' * wa), wa...))))
    qb = Quaternion(Δt / 2 * (qb * Quaternion(SVector(sqrt((2 / Δt)^2 - wb' * wb), wb...))))

    V12 * vrotate(SVector{3}(((xb + vrotate(SVector{3}(pb...), qb)) - (xa + vrotate(SVector{3}(pa...), qa)))...), inv(qa))
end

function transfunc1pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    v = vars[21:23]

    v' * vrotate(SVector{3}(((xb + vrotate(SVector{3}(pb...), qb)) - (xa + vrotate(SVector{3}(pa...), qa)))...), inv(qa))
end

function transfunc1vel(vars)
    Δt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    pa = vars[27:29]
    pb = vars[30:32]

    v = vars[33:35]

    xa = xa + Δt * va
    xb = xb + Δt * vb
    qa = Quaternion(Δt / 2 * (qa * Quaternion(SVector(sqrt((2 / Δt)^2 - wa' * wa), wa...))))
    qb = Quaternion(Δt / 2 * (qb * Quaternion(SVector(sqrt((2 / Δt)^2 - wb' * wb), wb...))))

    v' * vrotate(SVector{3}(((xb + vrotate(SVector{3}(pb...), qb)) - (xa + vrotate(SVector{3}(pa...), qa)))...), inv(qa))
end


function rotfunc3pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    offset = Quaternion(SVector(vars[15:18]...))

    
    VLᵀmat(offset) * Lᵀmat(qa) * qb
end

function rotfunc3vel(vars)
    Δt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    offset = Quaternion(SVector(vars[27:30]...))

    xa = xa + Δt * va
    xb = xb + Δt * vb
    qa = Quaternion(Δt / 2 * (qa * Quaternion(SVector(sqrt((2 / Δt)^2 - wa' * wa), wa...))))
    qb = Quaternion(Δt / 2 * (qb * Quaternion(SVector(sqrt((2 / Δt)^2 - wb' * wb), wb...))))

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

    V12 * (VLᵀmat(offset) * Lᵀmat(qa) * qb)
end

function rotfunc2vel(vars)
    Δt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    offset = Quaternion(SVector(vars[27:30]...))

    V1 = vars[31:33]
    V2 = vars[34:36]
    V12 = [V1';V2']

    xa = xa + Δt * va
    xb = xb + Δt * vb
    qa = Quaternion(Δt / 2 * (qa * Quaternion(SVector(sqrt((2 / Δt)^2 - wa' * wa), wa...))))
    qb = Quaternion(Δt / 2 * (qb * Quaternion(SVector(sqrt((2 / Δt)^2 - wb' * wb), wb...))))

    V12 * (VLᵀmat(offset) * Lᵀmat(qa) * qb)
end

function rotfunc1pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    offset = Quaternion(SVector(vars[15:18]...))

    V3 = vars[19:21]

    V3' * (VLᵀmat(offset) * Lᵀmat(qa) * qb)
end

function rotfunc1vel(vars)
    Δt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    offset = Quaternion(SVector(vars[27:30]...))

    V3 = vars[31:33]

    xa = xa + Δt * va
    xb = xb + Δt * vb
    qa = Quaternion(Δt / 2 * (qa * Quaternion(SVector(sqrt((2 / Δt)^2 - wa' * wa), wa...))))
    qb = Quaternion(Δt / 2 * (qb * Quaternion(SVector(sqrt((2 / Δt)^2 - wb' * wb), wb...))))

    V3' * (VLᵀmat(offset) * Lᵀmat(qa) * qb)
end

function transtest3()
    xa = rand(3)
    qa = Quaternion(rand(RotMatrix{3}))
    xb = rand(3)
    qb = Quaternion(rand(RotMatrix{3}))

    va = rand(3)
    wa = rand(3)
    vb = rand(3)
    wb = rand(3)

    pa = rand(3)
    pb = rand(3)


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational3(link1, link2, pa, pb))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])

    setPosition!(bot, link1, x = xa, q = qa)
    setPosition!(bot, link2, x = xb, q = qb)
    link1.s1 = SVector([va;wa]...)
    link2.s1 = SVector([vb;wb]...)

    res = ForwardDiff.jacobian(transfunc3pos, [xa;qa;xb;qb;pa;pb])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(joint1, 1, bot)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂pos(joint1, 1, bot)[1:3,4:6])
    n21 = norm(X2 - ∂g∂pos(joint1, 2, bot)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂pos(joint1, 2, bot)[1:3,4:6])

    res = ForwardDiff.jacobian(transfunc3vel, [xa;qa;xb;qb;va;wa;vb;wb;pa;pb])
    V1 = res[1:3,15:17]
    W1 = res[1:3,18:20]
    V2 = res[1:3,21:23]
    W2 = res[1:3,24:26]

    n31 = norm(V1 - ∂g∂vel(joint1, 1, bot)[1:3,1:3])
    n32 = norm(W1 - ∂g∂vel(joint1, 1, bot)[1:3,4:6])
    n41 = norm(V2 - ∂g∂vel(joint1, 2, bot)[1:3,1:3])
    n42 = norm(W2 - ∂g∂vel(joint1, 2, bot)[1:3,4:6])

    display((n11, n12, n21, n22))
    display((n31, n32, n41, n42))
    return
end


function transtest2()
    xa = rand(3)
    qa = Quaternion(rand(RotMatrix{3}))
    xb = rand(3)
    qb = Quaternion(rand(RotMatrix{3}))

    va = rand(3)
    wa = rand(3)
    vb = rand(3)
    wb = rand(3)

    pa = rand(3)
    pb = rand(3)

    v = rand(3)


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational2(link1, link2, pa, pb, v))
    V12 = joint1.constraints[1].V12

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])

    setPosition!(bot, link1, x = xa, q = qa)
    setPosition!(bot, link2, x = xb, q = qb)
    link1.s1 = SVector([va;wa]...)
    link2.s1 = SVector([vb;wb]...)


    res = ForwardDiff.jacobian(transfunc2pos, [xa;qa;xb;qb;pa;pb;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(joint1, 1, bot)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂pos(joint1, 1, bot)[1:2,4:6])
    n21 = norm(X2 - ∂g∂pos(joint1, 2, bot)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂pos(joint1, 2, bot)[1:2,4:6])


    res = ForwardDiff.jacobian(transfunc2vel, [xa;qa;xb;qb;va;wa;vb;wb;pa;pb;V12[1,:];V12[2,:]])
    V1 = res[1:2,15:17]
    W1 = res[1:2,18:20]
    V2 = res[1:2,21:23]
    W2 = res[1:2,24:26]

    n31 = norm(V1 - ∂g∂vel(joint1, 1, bot)[1:2,1:3])
    n32 = norm(W1 - ∂g∂vel(joint1, 1, bot)[1:2,4:6])
    n41 = norm(V2 - ∂g∂vel(joint1, 2, bot)[1:2,1:3])
    n42 = norm(W2 - ∂g∂vel(joint1, 2, bot)[1:2,4:6])

    display((n11, n12, n21, n22))
    display((n31, n32, n41, n42))
    return
end

function transtest1()
    xa = rand(3)
    qa = Quaternion(rand(RotMatrix{3}))
    xb = rand(3)
    qb = Quaternion(rand(RotMatrix{3}))

    va = rand(3)
    wa = rand(3)
    vb = rand(3)
    wb = rand(3)

    pa = rand(3)
    pb = rand(3)

    v = rand(3)


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational1(link1, link2, pa, pb, v))
    v = joint1.constraints[1].V3'

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])

    setPosition!(bot, link1, x = xa, q = qa)
    setPosition!(bot, link2, x = xb, q = qb)
    link1.s1 = SVector([va;wa]...)
    link2.s1 = SVector([vb;wb]...)


    res = ForwardDiff.gradient(transfunc1pos, [xa;qa;xb;qb;pa;pb;v])
    X1 = res[1:3]'
    Q1 = res[4:7]' * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[8:10]'
    Q2 = res[11:14]' * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(joint1, 1, bot)[1:3]')
    n12 = norm(Q1 - ∂g∂pos(joint1, 1, bot)[4:6]')
    n21 = norm(X2 - ∂g∂pos(joint1, 2, bot)[1:3]')
    n22 = norm(Q2 - ∂g∂pos(joint1, 2, bot)[4:6]')


    res = ForwardDiff.gradient(transfunc1vel, [xa;qa;xb;qb;va;wa;vb;wb;pa;pb;v])
    V1 = res[15:17]'
    W1 = res[18:20]'
    V2 = res[21:23]'
    W2 = res[24:26]'

    n31 = norm(V1 - ∂g∂vel(joint1, 1, bot)[1:3]')
    n32 = norm(W1 - ∂g∂vel(joint1, 1, bot)[4:6]')
    n41 = norm(V2 - ∂g∂vel(joint1, 2, bot)[1:3]')
    n42 = norm(W2 - ∂g∂vel(joint1, 2, bot)[4:6]')

    display((n11, n12, n21, n22))
    display((n31, n32, n41, n42))
    return
end


function rottest3()
    xa = rand(3)
    qa = Quaternion(rand(RotMatrix{3}))
    xb = rand(3)
    qb = Quaternion(rand(RotMatrix{3}))

    va = rand(3)
    wa = rand(3)
    vb = rand(3)
    wb = rand(3)

    offset = Quaternion(rand(RotMatrix{3}))


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational3(link1, link2, offset = offset))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])

    setPosition!(bot, link1, x = xa, q = qa)
    setPosition!(bot, link2, x = xb, q = qb)
    link1.s1 = SVector([va;wa]...)
    link2.s1 = SVector([vb;wb]...)


    res = ForwardDiff.jacobian(rotfunc3pos, [xa;qa;xb;qb;offset])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(joint1, 1, bot)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂pos(joint1, 1, bot)[1:3,4:6])
    n21 = norm(X2 - ∂g∂pos(joint1, 2, bot)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂pos(joint1, 2, bot)[1:3,4:6])


    res = ForwardDiff.jacobian(rotfunc3vel, [xa;qa;xb;qb;va;wa;vb;wb;offset])
    V1 = res[1:3,15:17]
    W1 = res[1:3,18:20]
    V2 = res[1:3,21:23]
    W2 = res[1:3,24:26]

    n31 = norm(V1 - ∂g∂vel(joint1, 1, bot)[1:3,1:3])
    n32 = norm(W1 - ∂g∂vel(joint1, 1, bot)[1:3,4:6])
    n41 = norm(V2 - ∂g∂vel(joint1, 2, bot)[1:3,1:3])
    n42 = norm(W2 - ∂g∂vel(joint1, 2, bot)[1:3,4:6])

    display((n11, n12, n21, n22))
    display((n31, n32, n41, n42))
    return
end

function rottest2()
    xa = rand(3)
    qa = Quaternion(rand(RotMatrix{3}))
    xb = rand(3)
    qb = Quaternion(rand(RotMatrix{3}))

    va = rand(3)
    wa = rand(3)
    vb = rand(3)
    wb = rand(3)

    offset = Quaternion(rand(RotMatrix{3}))

    v = rand(3)


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational2(link1, link2, v, offset = offset))
    V12 = joint1.constraints[1].V12

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])

    setPosition!(bot, link1, x = xa, q = qa)
    setPosition!(bot, link2, x = xb, q = qb)
    link1.s1 = SVector([va;wa]...)
    link2.s1 = SVector([vb;wb]...)


    res = ForwardDiff.jacobian(rotfunc2pos, [xa;qa;xb;qb;offset;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(joint1, 1, bot)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂pos(joint1, 1, bot)[1:2,4:6])
    n21 = norm(X2 - ∂g∂pos(joint1, 2, bot)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂pos(joint1, 2, bot)[1:2,4:6])


    res = ForwardDiff.jacobian(rotfunc2vel, [xa;qa;xb;qb;va;wa;vb;wb;offset;V12[1,:];V12[2,:]])
    V1 = res[1:2,15:17]
    W1 = res[1:2,18:20]
    V2 = res[1:2,21:23]
    W2 = res[1:2,24:26]

    n31 = norm(V1 - ∂g∂vel(joint1, 1, bot)[1:2,1:3])
    n32 = norm(W1 - ∂g∂vel(joint1, 1, bot)[1:2,4:6])
    n41 = norm(V2 - ∂g∂vel(joint1, 2, bot)[1:2,1:3])
    n42 = norm(W2 - ∂g∂vel(joint1, 2, bot)[1:2,4:6])

    display((n11, n12, n21, n22))
    display((n31, n32, n41, n42))
    return
end

function rottest1()
    xa = rand(3)
    qa = Quaternion(rand(RotMatrix{3}))
    xb = rand(3)
    qb = Quaternion(rand(RotMatrix{3}))

    va = rand(3)
    wa = rand(3)
    vb = rand(3)
    wb = rand(3)

    offset = Quaternion(rand(RotMatrix{3}))

    v = rand(3)


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational1(link1, link2, v, offset = offset))
    V3 = joint1.constraints[1].V3'

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])

    setPosition!(bot, link1, x = xa, q = qa)
    setPosition!(bot, link2, x = xb, q = qb)
    link1.s1 = SVector([va;wa]...)
    link2.s1 = SVector([vb;wb]...)


    res = ForwardDiff.gradient(rotfunc1pos, [xa;qa;xb;qb;offset;V3])
    X1 = res[1:3]'
    Q1 = res[4:7]' * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[8:10]'
    Q2 = res[11:14]' * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(joint1, 1, bot)[1:3]')
    n12 = norm(Q1 - ∂g∂pos(joint1, 1, bot)[4:6]')
    n21 = norm(X2 - ∂g∂pos(joint1, 2, bot)[1:3]')
    n22 = norm(Q2 - ∂g∂pos(joint1, 2, bot)[4:6]')


    res = ForwardDiff.gradient(rotfunc1vel, [xa;qa;xb;qb;va;wa;vb;wb;offset;V3])
    V1 = res[15:17]'
    W1 = res[18:20]'
    V2 = res[21:23]'
    W2 = res[24:26]'

    n31 = norm(V1 - ∂g∂vel(joint1, 1, bot)[1:3]')
    n32 = norm(W1 - ∂g∂vel(joint1, 1, bot)[4:6]')
    n41 = norm(V2 - ∂g∂vel(joint1, 2, bot)[1:3]')
    n42 = norm(W2 - ∂g∂vel(joint1, 2, bot)[4:6]')

    display((n11, n12, n21, n22))
    display((n31, n32, n41, n42))
    return
end