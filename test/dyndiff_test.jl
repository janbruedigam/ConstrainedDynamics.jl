using ConstrainedDynamics
using ForwardDiff
using Rotations
using StaticArrays
using LinearAlgebra

using ConstrainedDynamics: vrotate, Lmat, Vᵀmat, Lᵀmat, VLᵀmat, ∂g∂pos, ∂g∂vel, skew, skewplusdiag

function dynTvel(vars)
    ezg = SVector{3,T}(0, 0, 9.81)
    Δt = 0.01
    m = 1.

    v2 = vars[1:3]
    v1 = vars[4:6]

    m * ((v2 - v1) / Δt + ezg)
end

function dynTpos(vars)
    ezg = SVector{3,T}(0, 0, 9.81)
    Δt = 0.01
    m = 1.

    x3 = vars[1:3]
    x2 = vars[4:6]
    x1 = vars[7:9]

    v2 = (x3-x2)/Δt
    v1 = (x2-x1)/Δt

    m * ((v2 - v1) / Δt + ezg)
end

function dynRvel(vars)
    Δt = 0.01

    J = diagm(ones(3))
    ω2 = vars[1:3]
    ω1 = vars[4:6]
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    
    skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1)
end

function dynRpos(vars)
    Δt = 0.01

    J = diagm(ones(3))
    q3 = Quaternion(SVector(vars[1:4]...))
    q2 = Quaternion(SVector(vars[5:8]...))
    q1 = Quaternion(SVector(vars[9:12]...))

    ω2 = 2/Δt * VLᵀmat(q2)*q3
    ω1 = 2/Δt * VLᵀmat(q1)*q2
    sq2 = sqrt(4 / Δt^2 - ω2' * ω2)
    sq1 = sqrt(4 / Δt^2 - ω1' * ω1)
    
    skewplusdiag(ω2, sq2) * (J * ω2) - skewplusdiag(ω1, sq1) * (J * ω1)
end

function dyntestT()
    Δt = 0.01

    x3 = rand(3)
    x2 = rand(3)
    x1 = rand(3)
    q3 = Quaternion(rand(RotMatrix{3}))
    q2 = Quaternion(rand(RotMatrix{3}))
    q1 = Quaternion(rand(RotMatrix{3}))

    v2 = (x3-x2)/Δt
    v1 = (x2-x1)/Δt
    ω2 = 2/Δt * VLᵀmat(q2)*q3
    ω1 = 2/Δt * VLᵀmat(q1)*q2

    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link1.m = 1.0
    link1.J = diagm(ones(3))


    oc1 = EqualityConstraint(OriginConnection(origin, link1))

    mech = Mechanism(origin, [link1], [oc1])

    link1.x[1] = x1
    link1.x[2] = x2
    link1.q[1] = q1
    link1.q[2] = q2
    link1.s1 = SVector([v2;ω2]...)

    res = ForwardDiff.jacobian(dynTpos, [x3;x2;x1])
    X3 = res[1:3,1:3]

    n1 = norm(X3 - ∂dyn∂pos(link1, Δt)[1:3,1:3])

    res = ForwardDiff.jacobian(dynTvel, [v2;v1])
    V2 = res[1:3,1:3]

    n2 = norm(V2 - ∂dyn∂vel(link1, Δt)[1:3,1:3])

    # display(n1)
    # display(n2)
    return n1 + n2
end

function dyntestR()
    Δt = 0.01

    x3 = rand(3)
    x2 = rand(3)
    x1 = rand(3)
    q3 = Quaternion(rand(RotMatrix{3}))
    q2 = Quaternion(rand(RotMatrix{3}))
    q1 = Quaternion(rand(RotMatrix{3}))

    v2 = (x3-x2)/Δt
    v1 = (x2-x1)/Δt
    ω2 = 2/Δt * VLᵀmat(q2)*q3
    ω1 = 2/Δt * VLᵀmat(q1)*q2

    origin = Origin{Float64}()
    link1 = Body(Box(1., 1., 1., 1.))
    link1.m = 1.0
    link1.J = diagm(ones(3))


    oc1 = EqualityConstraint(OriginConnection(origin, link1))

    mech = Mechanism(origin, [link1], [oc1])

    link1.x[1] = x1
    link1.x[2] = x2
    link1.q[1] = q1
    link1.q[2] = q2
    link1.s1 = SVector([v2;ω2]...)

    res = ForwardDiff.jacobian(dynRpos, [q3;q2;q1])
    Q3 = res[1:3,1:4] * Lmat(Quaternion(q3)) * Vᵀmat()

    n1 = norm(Q3 - ∂dyn∂pos(link1, Δt)[4:6,4:6])

    res = ForwardDiff.jacobian(dynRvel, [ω2;ω1])
    W2 = res[1:3,1:3]

    n2 = norm(W2 - ∂dyn∂vel(link1, Δt)[4:6,4:6])

    # display(n1)
    # display(n2)
    return n1 + n2
end

for i=1:10
    @test isapprox(dyntestT(), 0.0; atol = 1e-8)
    @test isapprox(dyntestR(), 0.0; atol = 1e-8)
end