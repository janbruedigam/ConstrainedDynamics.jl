using ConstrainedDynamics
using Rotations
using StaticArrays

Δt = 0.01
length1 = 0.5
width, depth = 0.5, 0.5
box1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
box2 = Box(width, depth, length1, length1, color = RGBA(1., 0., 0.))

# Links
origin = Origin{Float64}()
link1 = Body(box1)
link2 = Body(box2)

# Constraints
joint1 = EqualityConstraint(OriginConnection(origin, link1))
joint2 = EqualityConstraint(OriginConnection(origin, link2))
# joint1 = EqualityConstraint(Spherical(origin, link1,zeros(3),zeros(3)))
# joint2 = EqualityConstraint(Spherical(link1, link2, [0;0;2.],zeros(3)))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box1;box2]


mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)

for i=1:10
    po1 = rand(3)
    p1o = rand(3)
    p12 = rand(3)
    p21 = rand(3)
    x1 = rand(3)
    Δx = rand(3)
    q1 = Quaternion(rand(RotMatrix{3}))
    Δq = Quaternion(rand(RotMatrix{3}))

    setPosition!(mech,link1, x = x1, q = q1, p1 = po1, p2 = p1o)
    setPosition!(mech,link1, link2, Δx = Δx, Δq = Δq, p1 = p12, p2 = p21)

    storage = simulate!(mech, 10., record = true)

    truex1 = [x1+po1+vrotate(SVector{3,Float64}(p1o),q1) for i=1:1000]
    trueq1 = [q1 for i=1:1000]
    truex2 = truex1 .+ [vrotate(SVector{3,Float64}(p12),q1) + Δx + vrotate(SVector{3,Float64}(p21),q2) for i=1:1000]
    trueq2 = [q1*Δq for i=1:1000]

    @test isapprox(sum(norm.(storage.x[1].-truex1)), 0.0; atol = 1e-8)
    @test isapprox(sum(norm.(storage.q[1].-trueq1)), 0.0; atol = 1e-8)
    @test isapprox(sum(norm.(storage.x[2].-truex2)), 0.0; atol = 1e-8)
    @test isapprox(sum(norm.(storage.q[2].-trueq2)), 0.0; atol = 1e-8)
end

for i=1:3
    axis1 = zeros(3)
    axis1[i] = 1
    for j=1:3
        axis2 = zeros(3)
        axis2[j] = 1

        po1 = rand(3)
        p1o = rand(3)
        p12 = rand(3)
        p21 = rand(3)
        x1 = rand(3)
        Δx = rand(3)
        q1 = Quaternion(rand(RotMatrix{3}))
        Δq = Quaternion(rand(RotMatrix{3}))
        v1 = rand(3)
        Δv = rand(3)

        setPosition!(mech,link1, x = x1, q = q1, p1 = po1, p2 = p1o)
        setPosition!(mech,link1, link2, Δx = Δx, Δq = Δq, p1 = p12, p2 = p21)
        
        setVelocity!(mech,link1,v = v1, ω = axis1, p1 = po1, p2 = p1o)
        setVelocity!(mech,link1,link2, Δv = Δv, Δω = axis2, p1 = p12, p2 = p21)

        storage = simulate!(mech, 10., record = true)

        truex10 = x1+po1+vrotate(SVector{3,Float64}(p1o),q1)
        trueq10 = q1
        truex20 = truex10 + vrotate(SVector{3,Float64}(p12),q1) + Δx + vrotate(SVector{3,Float64}(p21),q2)
        trueq20 = q1*Δq

        truex1 = [truex10 + v1*Δt*i for i=1:1000]
        trueq1 = [Quaternion(cos(i*Δt/2),(axis1*sin(i*Δt/2))...)*trueq10 for i=1:1000]
        # truex2 = [truex20 + (v1+Δv)*Δt*i for i=1:1000]
        # trueq2 = [Quaternion(cos(i*Δt/2),(axis1*sin(i*Δt/2))...)*trueq10 for i=1:1000]

        @test isapprox(sum(norm.(storage.x[1].-truex1)), 0.0; atol = 1e-8)
        @test isapprox(sum(norm.(storage.q[1].-trueq1)), 0.0; atol = 1e-8)
        # @test isapprox(sum(norm.(storage.x[2].-truex2)), 0.0; atol = 1e-8)
        # @test isapprox(sum(norm.(storage.q[2].-trueq2)), 0.0; atol = 1e-8)        
    end
end
