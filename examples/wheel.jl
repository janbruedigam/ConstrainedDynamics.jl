##########
# This is more of a test case than a proper working example
##########


using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]

l1 = 1.0 * 0.3
d = .05
b1 = Box(l1, d, d, l1, color = RGBA(1., 1., 0.))
b2 = Box(.1, 1., 1., 2., color = RGBA(1., 0., 0.))

vert11 = [l1 / 2;0.;0.]
vert12 = -vert11

# Initial orientation
phi1 = 0.
q1 = Quaternion(RotX(phi1))

# Links
origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b2)

# Constraints
socket0to1 = EqualityConstraint(Spherical(origin, link1, zeros(3), vert11))
joint1to5 = EqualityConstraint(Revolute(link1, link2, vert12, zeros(3), ex))

links = [link1;link2]
constraints = [socket0to1;joint1to5]
shapes = [b1;b2]


mech = Mechanism(origin, links, constraints;tend = 10.0,Δt = 0.001, shapes = shapes)
setPosition!(mech,origin,link1,p2 = vert11,Δq = q1)
setPosition!(mech,link1,link2,p1 = vert12,Δq = q1)
setVelocity!(mech,link1,ω = [50.;0;0])
setVelocity!(mech,link1,link2,p1 = vert12)

simulate!(mech,save = true)
visualize!(mech)
