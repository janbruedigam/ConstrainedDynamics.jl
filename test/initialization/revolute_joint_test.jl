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
    axis = rand(3)
    axis = axis/norm(axis)
    qoff = Quaternion{Float64}() # Quaternion(rand(RotMatrix{3}))


    # Constraints
    joint1 = EqualityConstraint(Revolute(origin, link1, p1, p2, axis, offset = qoff))

    links = [link1]
    constraints = [joint1]
    shapes = [box1]

    mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)

    xθ = rand(1)
    vω = rand(1)
    Fτ = rand(1)

    setPosition!(mech,joint1,SVector{1,Float64}(xθ))
    setVelocity!(mech,joint1,SVector{1,Float64}(vω))

    function control!(mechanism, k)
        if k==1
            setForce!(mech,joint1,SVector{1,Float64}(Fτ))
        else
            setForce!(mech,joint1,SVector{1,Float64}(Fτ)*0)
        end
        return
    end


    storage = simulate!(mech, 10., control!, record = true)

    angend = mod((xθ + (vω + Fτ*Δt)*10.0)[1],2pi)
    qend = qoff*Quaternion(cos(angend/2),(axis*sin(angend/2))...)
    minresult = minimalCoordinates(mech, joint1)
    if minresult[1] < 0
        minresult = minresult .+ 2pi
    end

    #TODO needs some velocity calculation fix, initial constraint satisfaction is not very good
    @test isapprox(norm(minresult - [angend]), 0.0; atol = 1.5*1e0)
    @test isapprox(norm(link1.state.xk[end] - (p1 - vrotate(SVector{3,Float64}(p2),qend))), 0.0; atol = 1.5*1e0)
    @test isapprox(norm((link1.state.qk[end]/qend)[2:4]), 0.0; atol = 1.5*1e0)
end

