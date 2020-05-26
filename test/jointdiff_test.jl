using ConstrainedDynamics
using ForwardDiff
using Rotations
using StaticArrays
using LinearAlgebra

using ConstrainedDynamics: vrotate, Lmat, Vᵀmat, Lᵀmat, VLᵀmat, ∂g∂pos, ∂g∂vel, skew, discretizestate!, 
    transfunc1pos, transfunc2pos, transfunc3pos, transfunc1vel, transfunc2vel, transfunc3vel, 
    rotfunc1pos, rotfunc2pos, rotfunc3pos, rotfunc1vel, rotfunc2vel, rotfunc3vel, 
    getx1, getq1, getx2, getq2, getv1, getω1, getv2, getω2, getxd3, getqd3, getxd2, getqd2


function transtest3()
    Δt = 0.01
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
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational3{Float64}(body1, body2, p1=pa, p2=pb))

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)
    body1.state.vc[2] = va
    body1.state.ωc[2] = wa
    body2.state.vc[2] = vb
    body2.state.ωc[2] = wb

    res = ForwardDiff.jacobian(transfunc3pos, [getxd2(body1);getqd2(body1);getxd2(body2);getqd2(body2);pa;pb])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * Lmat(Quaternion(getqd2(body1))) * Vᵀmat()
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * Lmat(Quaternion(getqd2(body2))) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:3,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:3,4:6])

    res = ForwardDiff.jacobian(transfunc3vel, [getxd2(body1);getqd2(body1);getv2(body1);getω2(body1);getxd2(body2);getqd2(body2);getv2(body2);getω2(body2);pa;pb])
    V1 = res[1:3,8:10]
    W1 = res[1:3,11:13]
    V2 = res[1:3,21:23]
    W2 = res[1:3,24:26]

    n31 = norm(V1 - ∂g∂vel(mech, joint1, 1)[1:3,1:3])
    n32 = norm(W1 - ∂g∂vel(mech, joint1, 1)[1:3,4:6])
    n41 = norm(V2 - ∂g∂vel(mech, joint1, 2)[1:3,1:3])
    n42 = norm(W2 - ∂g∂vel(mech, joint1, 2)[1:3,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end


function transtest2()
    Δt = 0.01
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
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational2{Float64}(body1, body2, p1=pa, p2=pb, axis=v))
    V12 = joint1.constraints[1].V12

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)
    body1.state.vc[2] = va
    body1.state.ωc[2] = wa
    body2.state.vc[2] = vb
    body2.state.ωc[2] = wb

    res = ForwardDiff.jacobian(transfunc2pos, [getxd2(body1);getqd2(body1);getxd2(body2);getqd2(body2);pa;pb;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * Lmat(Quaternion(getqd2(body1))) * Vᵀmat()
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * Lmat(Quaternion(getqd2(body2))) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:2,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:2,4:6])


    res = ForwardDiff.jacobian(transfunc2vel, [getxd2(body1);getqd2(body1);getv2(body1);getω2(body1);getxd2(body2);getqd2(body2);getv2(body2);getω2(body2);pa;pb;V12[1,:];V12[2,:]])
    V1 = res[1:2,8:10]
    W1 = res[1:2,11:13]
    V2 = res[1:2,21:23]
    W2 = res[1:2,24:26]

    n31 = norm(V1 - ∂g∂vel(mech, joint1, 1)[1:2,1:3])
    n32 = norm(W1 - ∂g∂vel(mech, joint1, 1)[1:2,4:6])
    n41 = norm(V2 - ∂g∂vel(mech, joint1, 2)[1:2,1:3])
    n42 = norm(W2 - ∂g∂vel(mech, joint1, 2)[1:2,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

function transtest1()
    Δt = 0.01
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
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational1{Float64}(body1, body2, p1=pa, p2=pb, axis=v))
    v = joint1.constraints[1].V3'

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)
    body1.state.vc[2] = va
    body1.state.ωc[2] = wa
    body2.state.vc[2] = vb
    body2.state.ωc[2] = wb

    res = ForwardDiff.gradient(transfunc1pos, [getxd2(body1);getqd2(body1);getxd2(body2);getqd2(body2);pa;pb;v])
    X1 = res[1:3]'
    Q1 = res[4:7]' * Lmat(Quaternion(getqd2(body1))) * Vᵀmat()
    X2 = res[8:10]'
    Q2 = res[11:14]' * Lmat(Quaternion(getqd2(body2))) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3]')
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[4:6]')
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3]')
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[4:6]')


    res = ForwardDiff.gradient(transfunc1vel, [getxd2(body1);getqd2(body1);getv2(body1);getω2(body1);getxd2(body2);getqd2(body2);getv2(body2);getω2(body2);pa;pb;v])
    V1 = res[8:10]'
    W1 = res[11:13]'
    V2 = res[21:23]'
    W2 = res[24:26]'

    n31 = norm(V1 - ∂g∂vel(mech, joint1, 1)[1:3]')
    n32 = norm(W1 - ∂g∂vel(mech, joint1, 1)[4:6]')
    n41 = norm(V2 - ∂g∂vel(mech, joint1, 2)[1:3]')
    n42 = norm(W2 - ∂g∂vel(mech, joint1, 2)[4:6]')

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end


function rottest3()
    Δt = 0.01
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
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational3{Float64}(body1, body2, offset = offset))

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)
    body1.state.vc[2] = va
    body1.state.ωc[2] = wa
    body2.state.vc[2] = vb
    body2.state.ωc[2] = wb

    res = ForwardDiff.jacobian(rotfunc3pos, [getxd2(body1);getqd2(body1);getxd2(body2);getqd2(body2);offset])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:3,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:3,4:6])


    res = ForwardDiff.jacobian(rotfunc3vel, [getxd2(body1);getqd2(body1);getv2(body1);getω2(body1);getxd2(body2);getqd2(body2);getv2(body2);getω2(body2);offset])
    V1 = res[1:3,8:10]
    W1 = res[1:3,11:13]
    V2 = res[1:3,21:23]
    W2 = res[1:3,24:26]

    n31 = norm(V1 - ∂g∂vel(mech, joint1, 1)[1:3,1:3])
    n32 = norm(W1 - ∂g∂vel(mech, joint1, 1)[1:3,4:6])
    n41 = norm(V2 - ∂g∂vel(mech, joint1, 2)[1:3,1:3])
    n42 = norm(W2 - ∂g∂vel(mech, joint1, 2)[1:3,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

function rottest2()
    Δt = 0.01
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
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational2{Float64}(body1, body2, axis = v, offset = offset))
    V12 = joint1.constraints[1].V12

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)
    body1.state.vc[2] = va
    body1.state.ωc[2] = wa
    body2.state.vc[2] = vb
    body2.state.ωc[2] = wb

    res = ForwardDiff.jacobian(rotfunc2pos, [getxd2(body1);getqd2(body1);getxd2(body2);getqd2(body2);offset;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:2,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:2,4:6])


    res = ForwardDiff.jacobian(rotfunc2vel, [getxd2(body1);getqd2(body1);getv2(body1);getω2(body1);getxd2(body2);getqd2(body2);getv2(body2);getω2(body2);offset;V12[1,:];V12[2,:]])
    V1 = res[1:2,8:10]
    W1 = res[1:2,11:13]
    V2 = res[1:2,21:23]
    W2 = res[1:2,24:26]

    n31 = norm(V1 - ∂g∂vel(mech, joint1, 1)[1:2,1:3])
    n32 = norm(W1 - ∂g∂vel(mech, joint1, 1)[1:2,4:6])
    n41 = norm(V2 - ∂g∂vel(mech, joint1, 2)[1:2,1:3])
    n42 = norm(W2 - ∂g∂vel(mech, joint1, 2)[1:2,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

function rottest1()
    Δt = 0.01
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
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational1{Float64}(body1, body2, axis = v, offset = offset))
    V3 = joint1.constraints[1].V3'

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)
    body1.state.vc[2] = va
    body1.state.ωc[2] = wa
    body2.state.vc[2] = vb
    body2.state.ωc[2] = wb

    res = ForwardDiff.gradient(rotfunc1pos, [getxd2(body1);getqd2(body1);getxd2(body2);getqd2(body2);offset;V3])
    X1 = res[1:3]'
    Q1 = res[4:7]' * Lmat(Quaternion(qa)) * Vᵀmat()
    X2 = res[8:10]'
    Q2 = res[11:14]' * Lmat(Quaternion(qb)) * Vᵀmat()

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3]')
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[4:6]')
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3]')
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[4:6]')


    res = ForwardDiff.gradient(rotfunc1vel, [getxd2(body1);getqd2(body1);getv2(body1);getω2(body1);getxd2(body2);getqd2(body2);getv2(body2);getω2(body2);offset;V3])
    V1 = res[8:10]'
    W1 = res[11:13]'
    V2 = res[21:23]'
    W2 = res[24:26]'

    n31 = norm(V1 - ∂g∂vel(mech, joint1, 1)[1:3]')
    n32 = norm(W1 - ∂g∂vel(mech, joint1, 1)[4:6]')
    n41 = norm(V2 - ∂g∂vel(mech, joint1, 2)[1:3]')
    n42 = norm(W2 - ∂g∂vel(mech, joint1, 2)[4:6]')

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

for i=1:10
    @test isapprox(transtest3(), 0.0; atol = 1e-8)
    @test isapprox(transtest2(), 0.0; atol = 1e-8)
    @test isapprox(transtest1(), 0.0; atol = 1e-8)
    @test isapprox(rottest3(), 0.0; atol = 1e-8)
    @test isapprox(rottest2(), 0.0; atol = 1e-8)
    @test isapprox(rottest1(), 0.0; atol = 1e-8)
end