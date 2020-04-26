using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 1.0
x, z = .1, .1
b1 = Box(x, l1, z, l1 * 100, color = RGBA(1., 1., 0.))
b2 = Box(x, l1, z, l1 * 1, color = RGBA(1., 0., 0.))

vert11 = [0.;l1 / 2;0.]
vert12 = -vert11


# Links
origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b2)

# Constraints
joint0to1 = EqualityConstraint(OriginConnection(origin, link1))
# joint0to1 = EqualityConstraint(Fixed(origin,link1,zeros(3),zeros(3)))
joint1to2 = EqualityConstraint(Revolute(link1, link2, vert12, vert11, ex))


links = [link1; link2]
constraints = [joint0to1;joint1to2]
shapes = [b1;b2]

mech = Mechanism(origin, links, constraints, shapes = shapes, g = 0.)
setPosition!(mech,link1,link2,p1 = vert12,p2 = vert11)

setForce!(mech, joint1to2, [nothing;[0.05]])

simulate!(mech,save = true)
visualize!(mech)
