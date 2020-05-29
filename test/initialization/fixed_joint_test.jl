using ConstrainedDynamics
using ConstrainedDynamics: vrotate
using Rotations
using StaticArrays
using LinearAlgebra


Δt = 0.01
length1 = 1.0
width, depth = 1.0, 1.0
box1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

for i=1:10
    # Links
    origin = Origin{Float64}()
    link1 = Body(box1)
    link1.m = 1.0
    link1.J = diagm(ones(3))

    p1 = rand(3)
    p2 = rand(3)
    qoff = Quaternion(rand(RotMatrix{3}))


    # Constraints
    joint1 = EqualityConstraint(Fixed(origin, link1, p1, p2, offset = qoff))

    links = [link1]
    constraints = [joint1]
    shapes = [box1]

    mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)

    setPosition!(mech,joint1,SVector{0,Float64}())
    setVelocity!(mech,joint1,SVector{0,Float64}())
    setForce!(mech,joint1,SVector{0,Float64}())

    storage = simulate!(mech, 10., record = true)

    @test isapprox(norm(minimalCoordinates(mech, joint1) - SVector{0,Float64}()), 0.0; atol = 1e-8)
    @test isapprox(norm(link1.state.xk[end] - (p1 - vrotate(SVector{3,Float64}(p2),qoff))), 0.0; atol = 1e-8)
    @test isapprox(norm(link1.state.qk[end] - qoff), 0.0; atol = 1e-8)
end