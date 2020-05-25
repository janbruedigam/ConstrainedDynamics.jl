using ConstrainedDynamics
using ForwardDiff
using Rotations
using StaticArrays
using LinearAlgebra

using ConstrainedDynamics: vrotate, Lmat, Vᵀmat, Lᵀmat, VLᵀmat, ∂dyn∂pos, ∂dyn∂vel, skew, skewplusdiag, discretizestate!,
    velTpos, dynTvel, velRpos, dynRvel


function dyntestT()
    Δt = 0.01

    x1 = rand(3)
    v1 = rand(3)
    v2 = rand(3)

    q1 = Quaternion(rand(RotMatrix{3}))
    ω1 = @SVector rand(3)
    ω2 = @SVector rand(3)


    origin = Origin{Float64}()
    body1 = Body(Box(1., 1., 1., 1.))
    body1.m = 1.0
    body1.J = diagm(ones(3))


    oc1 = EqualityConstraint(OriginConnection(origin, body1))

    mech = Mechanism(origin, [body1], [oc1])
    discretizestate!(body1,x1,q1,v1,v2,ω1,ω2,Δt)


    res = ForwardDiff.jacobian(dynTvel, [v2;v1])
    X2 = res[1:3,1:3] * velTpos([x1;v2])

    n1 = norm(X2 - ∂dyn∂pos(body1, Δt)[1:3,1:3])

    res = ForwardDiff.jacobian(dynTvel, [v2;v1])
    V2 = res[1:3,1:3]

    n2 = norm(V2 - ∂dyn∂vel(body1, Δt)[1:3,1:3])

    # display(n1)
    # display(n2)
    return n1 + n2
end

function dyntestR()
    Δt = 0.01

    x1 = rand(3)
    v1 = rand(3)
    v2 = rand(3)

    q1 = Quaternion(rand(RotMatrix{3}))
    ω1 = @SVector rand(3)
    ω2 = @SVector rand(3)


    origin = Origin{Float64}()
    body1 = Body(Box(1., 1., 1., 1.))
    body1.m = 1.0
    body1.J = diagm(ones(3))


    oc1 = EqualityConstraint(OriginConnection(origin, body1))

    mech = Mechanism(origin, [body1], [oc1])
    discretizestate!(body1,x1,q1,v1,v2,ω1,ω2,Δt)

    res = ForwardDiff.jacobian(dynRvel, [ω2;ω1])
    Q2 = res[1:3,1:3] * velRpos([q1;ω2])

    n1 = norm(Q2 - ∂dyn∂pos(body1, Δt)[4:6,4:6])

    res = ForwardDiff.jacobian(dynRvel, [ω2;ω1])
    W2 = res[1:3,1:3]

    n2 = norm(W2 - ∂dyn∂vel(body1, Δt)[4:6,4:6])

    # display(n1)
    # display(n2)
    return n1 + n2
end

for i=1:10
    @test isapprox(dyntestT(), 0.0; atol = 1e-8)
    @test isapprox(dyntestR(), 0.0; atol = 1e-8)
end