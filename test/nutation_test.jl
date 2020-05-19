using ConstrainedDynamics
using LinearAlgebra

# Parameters
h = .1
r = 1.
b1 = Cylinder(r, h, h*r, color = RGBA(1., 0., 0.))

# length1 = 0.1
# width, depth = 2., 2.
# b1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

# Links
origin = Origin{Float64}()
link1 = Body(b1)

# Constraints
joint1 = EqualityConstraint(OriginConnection(origin, link1))

links = [link1]
constraints = [joint1]
shapes = [b1]


mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)

axis = [0;0;1.]
setVelocity!(mech,link1, ω = 2pi*axis)

function control!(mechanism, k)
    if k==1
        setForce!(mechanism, mechanism.bodies[1], F = [0;0;2.], r=[0;1.;0])
    else
        setForce!(mechanism, mechanism.bodies[1], F = [0;0;0.], r=[0;0.0;0])
    end
    return
end

storage = simulate!(mech, 10., controller!, record = true)