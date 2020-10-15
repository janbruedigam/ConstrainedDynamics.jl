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
joint1 = EqualityConstraint(Floating(origin, link1))

links = [link1]
constraints = [joint1]
shapes = [b1]


mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)

axis = [0;0;1.]
speed = 2pi #*0
setVelocity!(link1, ω = speed*axis)

function nutation_control!(mechanism, k)
    if k==1
        setForce!(mechanism.bodies[1], F = SA[0;0;0.2] * 0, p = SA[0;1.;0])
    else
        setForce!(mechanism.bodies[1], F = szeros(3), p = szeros(3))
    end
    return
end
