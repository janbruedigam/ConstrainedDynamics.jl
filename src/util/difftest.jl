using ForwardDiff
using Rotations
using StaticArrays
using LinearAlgebra

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics
using Main.MaximalCoordinateDynamics: vrotate, Lmat, Vᵀmat, ∂g∂pos, ∂g∂vel, skew

function transfunc0pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    (xb + vrotate(pb, qb)) - (xa + vrotate(pa, qa))
end

function transfunc0vel(vars)
    dt = 0.01

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

    xa = xa + dt * va
    xb = xb + dt * vb
    qa = Quaternion(dt / 2 * (qa * Quaternion(SVector(sqrt((2 / dt)^2 - wa' * wa), wa...))))
    qb = Quaternion(dt / 2 * (qb * Quaternion(SVector(sqrt((2 / dt)^2 - wb' * wb), wb...))))

    (xb + vrotate(pb, qb)) - (xa + vrotate(pa, qa))
end

function transfunc1pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    V1 = vars[21:23]
    V2 = vars[24:26]
    V12 = [V1';V2']

    V12 * vrotate((xb + vrotate(pb, qb)) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc1vel(vars)
    dt = 0.01

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

    xa = xa + dt * va
    xb = xb + dt * vb
    qa = Quaternion(dt / 2 * (qa * Quaternion(SVector(sqrt((2 / dt)^2 - wa' * wa), wa...))))
    qb = Quaternion(dt / 2 * (qb * Quaternion(SVector(sqrt((2 / dt)^2 - wb' * wb), wb...))))

    V12 * vrotate((xb + vrotate(pb, qb)) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc2pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    pa = vars[15:17]
    pb = vars[18:20]

    v = vars[21:23]

    v' * vrotate((xb + vrotate(pb, qb)) - (xa + vrotate(pa, qa)), inv(qa))
end

function transfunc2vel(vars)
    dt = 0.01

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

    xa = xa + dt * va
    xb = xb + dt * vb
    qa = Quaternion(dt / 2 * (qa * Quaternion(SVector(sqrt((2 / dt)^2 - wa' * wa), wa...))))
    qb = Quaternion(dt / 2 * (qb * Quaternion(SVector(sqrt((2 / dt)^2 - wb' * wb), wb...))))

    v' * vrotate((xb + vrotate(pb, qb)) - (xa + vrotate(pa, qa)), inv(qa))
end


function rotfunc0pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    offset = Quaternion(SVector(vars[15:18]...))

    
    MaximalCoordinateDynamics.VLᵀmat(qa) * qb - offset[2:4]
end

function rotfunc0vel(vars)
    dt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    offset = Quaternion(SVector(vars[27:30]...))

    xa = xa + dt * va
    xb = xb + dt * vb
    qa = Quaternion(dt / 2 * (qa * Quaternion(SVector(sqrt((2 / dt)^2 - wa' * wa), wa...))))
    qb = Quaternion(dt / 2 * (qb * Quaternion(SVector(sqrt((2 / dt)^2 - wb' * wb), wb...))))

    MaximalCoordinateDynamics.VLᵀmat(qa) * qb - offset[2:4]
end

function rotfunc1pos(vars)
    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    V1 = vars[15:17]
    V2 = vars[18:20]
    V12 = [V1';V2']

    V12 * (MaximalCoordinateDynamics.VLᵀmat(qa) * qb)
end

function rotfunc1vel(vars)
    dt = 0.01

    xa = vars[1:3]
    qa = Quaternion(SVector(vars[4:7]...))
    xb = vars[8:10]
    qb = Quaternion(SVector(vars[11:14]...))

    va = vars[15:17]
    wa = vars[18:20]
    vb = vars[21:23]
    wb = vars[24:26]

    V1 = vars[27:29]
    V2 = vars[30:32]
    V12 = [V1';V2']

    xa = xa + dt * va
    xb = xb + dt * vb
    qa = Quaternion(dt / 2 * (qa * Quaternion(SVector(sqrt((2 / dt)^2 - wa' * wa), wa...))))
    qb = Quaternion(dt / 2 * (qb * Quaternion(SVector(sqrt((2 / dt)^2 - wb' * wb), wb...))))

    V12 * (MaximalCoordinateDynamics.VLᵀmat(qa) * qb)
end

function transtest0()
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
    setInit!(origin, link1, xa, zeros(3), q = qa)
    link1.s1 = SVector([va;wa]...)
    link2 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link2, xb, zeros(3), q = qb)
    link2.s1 = SVector([vb;wb]...)

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(MaximalCoordinateDynamics.Translational0(link1, link2, pa, pb))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])


    res = ForwardDiff.jacobian(transfunc0pos, [xa;qa;xb;qb;pa;pb])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n1 = norm([X1 Q1] - ∂g∂pos(joint1, 1, bot))
    n2 = norm([X2 Q2] - ∂g∂pos(joint1, 2, bot))


    res = ForwardDiff.jacobian(transfunc0vel, [xa;qa;xb;qb;va;wa;vb;wb;pa;pb])
    V1 = res[1:3,15:17]
    W1 = res[1:3,18:20]
    V2 = res[1:3,21:23]
    W2 = res[1:3,24:26]

    n3 = norm([V1 W1] - ∂g∂vel(joint1, 1, bot))
    n4 = norm([V2 W2] - ∂g∂vel(joint1, 2, bot))

    return n1, n2, n3, n4
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
    V12 = (@SMatrix [1 0 0; 0 1 0]) * svd(skew(v)).Vt


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link1, xa, zeros(3), q = qa)
    link1.s1 = SVector([va;wa]...)
    link2 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link2, xb, zeros(3), q = qb)
    link2.s1 = SVector([vb;wb]...)

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(MaximalCoordinateDynamics.Translational1(link1, link2, pa, pb, v))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])


    res = ForwardDiff.jacobian(transfunc1pos, [xa;qa;xb;qb;pa;pb;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n1 = norm([X1 Q1] - ∂g∂pos(joint1, 1, bot))
    n2 = norm([X2 Q2] - ∂g∂pos(joint1, 2, bot))


    res = ForwardDiff.jacobian(transfunc1vel, [xa;qa;xb;qb;va;wa;vb;wb;pa;pb;V12[1,:];V12[2,:]])
    V1 = res[1:2,15:17]
    W1 = res[1:2,18:20]
    V2 = res[1:2,21:23]
    W2 = res[1:2,24:26]

    n3 = norm([V1 W1] - ∂g∂vel(joint1, 1, bot))
    n4 = norm([V2 W2] - ∂g∂vel(joint1, 2, bot))

    return n1, n2, n3, n4
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
    setInit!(origin, link1, xa, zeros(3), q = qa)
    link1.s1 = SVector([va;wa]...)
    link2 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link2, xb, zeros(3), q = qb)
    link2.s1 = SVector([vb;wb]...)

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(MaximalCoordinateDynamics.Translational2(link1, link2, pa, pb, v))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])


    res = ForwardDiff.gradient(transfunc2pos, [xa;qa;xb;qb;pa;pb;v])
    X1 = res[1:3]'
    Q1 = res[4:7]' * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[8:10]'
    Q2 = res[11:14]' * Lmat(Quaternion(qb)) * Vᵀmat()

    n1 = norm([X1 Q1] - ∂g∂pos(joint1, 1, bot))
    n2 = norm([X2 Q2] - ∂g∂pos(joint1, 2, bot))


    res = ForwardDiff.gradient(transfunc2vel, [xa;qa;xb;qb;va;wa;vb;wb;pa;pb;v])
    V1 = res[15:17]'
    W1 = res[18:20]'
    V2 = res[21:23]'
    W2 = res[24:26]'

    n3 = norm([V1 W1] - ∂g∂vel(joint1, 1, bot))
    n4 = norm([V2 W2] - ∂g∂vel(joint1, 2, bot))

    return n1, n2, n3, n4
end


function rottest0()
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
    setInit!(origin, link1, xa, zeros(3), q = qa)
    link1.s1 = SVector([va;wa]...)
    link2 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link2, xb, zeros(3), q = qb)
    link2.s1 = SVector([vb;wb]...)

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(MaximalCoordinateDynamics.Rotational0(link1, link2, offset=offset))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])


    res = ForwardDiff.jacobian(rotfunc0pos, [xa;qa;xb;qb;offset])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n1 = norm([X1 Q1] - ∂g∂pos(joint1, 1, bot))
    n2 = norm([X2 Q2] - ∂g∂pos(joint1, 2, bot))


    res = ForwardDiff.jacobian(rotfunc0vel, [xa;qa;xb;qb;va;wa;vb;wb;offset])
    V1 = res[1:3,15:17]
    W1 = res[1:3,18:20]
    V2 = res[1:3,21:23]
    W2 = res[1:3,24:26]

    n3 = norm([V1 W1] - ∂g∂vel(joint1, 1, bot))
    n4 = norm([V2 W2] - ∂g∂vel(joint1, 2, bot))

    return n1, n2, n3, n4
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

    v = rand(3)
    V12 = (@SMatrix [1 0 0; 0 1 0]) * svd(skew(v)).Vt


    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link1, xa, zeros(3), q = qa)
    link1.s1 = SVector([va;wa]...)
    link2 = Body(Box(1., 1., 1., 1.))
    setInit!(origin, link2, xb, zeros(3), q = qb)
    link2.s1 = SVector([vb;wb]...)

    oc1 = EqualityConstraint(OriginConnection(origin, link1))
    oc2 = EqualityConstraint(OriginConnection(origin, link2))
    joint1 = EqualityConstraint(MaximalCoordinateDynamics.Rotational1(link1, link2, v))

    bot = Mechanism(origin, [link1;link2], [oc1;oc2;joint1])


    res = ForwardDiff.jacobian(rotfunc1pos, [xa;qa;xb;qb;;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n1 = norm([X1 Q1] - ∂g∂pos(joint1, 1, bot))
    n2 = norm([X2 Q2] - ∂g∂pos(joint1, 2, bot))


    res = ForwardDiff.jacobian(rotfunc1vel, [xa;qa;xb;qb;va;wa;vb;wb;;V12[1,:];V12[2,:]])
    V1 = res[1:2,15:17]
    W1 = res[1:2,18:20]
    V2 = res[1:2,21:23]
    W2 = res[1:2,24:26]

    n3 = norm([V1 W1] - ∂g∂vel(joint1, 1, bot))
    n4 = norm([V2 W2] - ∂g∂vel(joint1, 2, bot))

    return n1, n2, n3, n4
end