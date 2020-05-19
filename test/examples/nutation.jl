using ConstrainedDynamics
using LinearAlgebra
using StaticArrays


# Parameters
h = .1
r = 1.
b1 = Cylinder(r, h, h*r, color = RGBA(1., 0., 0.))

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
speed = 20pi #*0
setVelocity!(mech,link1, ω = speed*axis)

function control!(mechanism, k)
    if k==1
        setForce!(mechanism, mechanism.bodies[1], F = SA[0;0;2.], r=SA[0;1.;0])
    else
        setForce!(mechanism, mechanism.bodies[1], F = SA[0;0;0.], r=SA[0;0.0;0])
    end
    return
end
