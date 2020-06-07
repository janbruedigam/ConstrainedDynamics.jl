using ConstrainedDynamics
using ForwardDiff
using Rotations
using Rotations: params
using StaticArrays
using LinearAlgebra

using ConstrainedDynamics: ∂g∂pos, ∂g∂vel, discretizestate!, 
    transfunc1pos, transfunc2pos, transfunc3pos, transfunc1vel, transfunc2vel, transfunc3vel, 
    rotfunc1pos, rotfunc2pos, rotfunc3pos, rotfunc1vel, rotfunc2vel, rotfunc3vel, 
    getxqkvector, getxk, getqk, getstateandvestimate, LVᵀmat


function transtest3()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = Quaternion(rand(RotMatrix{3}))
    xb1 = rand(3)
    qb1 = Quaternion(rand(RotMatrix{3}))

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    pa = rand(3)
    pb = rand(3)


    origin = Origin{Float64}()
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational3{Float64}(body1, body2, p1=pa, p2=pb))

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.jacobian(transfunc3pos, [getxqkvector(body1.state);getxqkvector(body2.state);pa;pb])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:3,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:3,4:6])

    res = ForwardDiff.jacobian(transfunc3vel, [getstateandvestimate(body1.state);getstateandvestimate(body2.state);pa;pb])
    V1 = res[1:3,14:16]
    W1 = res[1:3,17:19]
    V2 = res[1:3,33:35]
    W2 = res[1:3,36:38]

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

    res = ForwardDiff.jacobian(transfunc2pos, [getxqkvector(body1.state);getxqkvector(body2.state);pa;pb;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:2,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:2,4:6])


    res = ForwardDiff.jacobian(transfunc2vel, [getstateandvestimate(body1.state);getstateandvestimate(body2.state);pa;pb;V12[1,:];V12[2,:]])
    V1 = res[1:2,14:16]
    W1 = res[1:2,17:19]
    V2 = res[1:2,33:35]
    W2 = res[1:2,36:38]

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

    res = ForwardDiff.gradient(transfunc1pos, [getxqkvector(body1.state);getxqkvector(body2.state);pa;pb;v])
    X1 = res[1:3]'
    Q1 = res[4:7]' * LVᵀmat(getqk(body1.state))
    X2 = res[8:10]'
    Q2 = res[11:14]' * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3]')
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[4:6]')
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3]')
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[4:6]')


    res = ForwardDiff.gradient(transfunc1vel, [getstateandvestimate(body1.state);getstateandvestimate(body2.state);pa;pb;v])
    V1 = res[14:16]'
    W1 = res[17:19]'
    V2 = res[33:35]'
    W2 = res[36:38]'

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

    qoffset = Quaternion(rand(RotMatrix{3}))


    origin = Origin{Float64}()
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational3{Float64}(body1, body2, qoffset = qoffset))

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)

    res = ForwardDiff.jacobian(rotfunc3pos, [getxqkvector(body1.state);getxqkvector(body2.state);params(qoffset)])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:3,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:3,4:6])


    res = ForwardDiff.jacobian(rotfunc3vel, [getstateandvestimate(body1.state);getstateandvestimate(body2.state);params(qoffset)])
    V1 = res[1:3,14:16]
    W1 = res[1:3,17:19]
    V2 = res[1:3,33:35]
    W2 = res[1:3,36:38]

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

    qoffset = Quaternion(rand(RotMatrix{3}))

    v = rand(3)


    origin = Origin{Float64}()
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational2{Float64}(body1, body2, axis = v, qoffset = qoffset))
    V12 = joint1.constraints[1].V12

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)

    res = ForwardDiff.jacobian(rotfunc2pos, [getxqkvector(body1.state);getxqkvector(body2.state);params(qoffset);V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[1:2,4:6])
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[1:2,4:6])


    res = ForwardDiff.jacobian(rotfunc2vel, [getstateandvestimate(body1.state);getstateandvestimate(body2.state);params(qoffset);V12[1,:];V12[2,:]])
    V1 = res[1:2,14:16]
    W1 = res[1:2,17:19]
    V2 = res[1:2,33:35]
    W2 = res[1:2,36:38]

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

    qoffset = Quaternion(rand(RotMatrix{3}))

    v = rand(3)


    origin = Origin{Float64}()
    body1 = Body(Box(1., 1., 1., 1.))
    body2 = Body(Box(1., 1., 1., 1.))

    oc1 = EqualityConstraint(OriginConnection(origin, body1))
    oc2 = EqualityConstraint(OriginConnection(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational1{Float64}(body1, body2, axis = v, qoffset = qoffset))
    V3 = joint1.constraints[1].V3'

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])

    setPosition!(body1, x = xa, q = qa)
    setPosition!(body2, x = xb, q = qb)
    discretizestate!(body1, Δt)
    discretizestate!(body2, Δt)

    res = ForwardDiff.gradient(rotfunc1pos, [getxqkvector(body1.state);getxqkvector(body2.state);params(qoffset);V3])
    X1 = res[1:3]'
    Q1 = res[4:7]' * LVᵀmat(getqk(body1.state))
    X2 = res[8:10]'
    Q2 = res[11:14]' * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂pos(mech, joint1, 1)[1:3]')
    n12 = norm(Q1 - ∂g∂pos(mech, joint1, 1)[4:6]')
    n21 = norm(X2 - ∂g∂pos(mech, joint1, 2)[1:3]')
    n22 = norm(Q2 - ∂g∂pos(mech, joint1, 2)[4:6]')


    res = ForwardDiff.gradient(rotfunc1vel, [getstateandvestimate(body1.state);getstateandvestimate(body2.state);params(qoffset);V3])
    V1 = res[14:16]'
    W1 = res[17:19]'
    V2 = res[33:35]'
    W2 = res[36:38]'

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
