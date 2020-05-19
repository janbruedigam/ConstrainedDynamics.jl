using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5
box1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
box2 = Box(width, depth, length1, length1, color = RGBA(1., 0., 0.))

# Links
origin = Origin{Float64}()
link1 = Body(box1)
link2 = Body(box2)

# Constraints
# joint1 = EqualityConstraint(OriginConnection(origin, link1))
# joint2 = EqualityConstraint(OriginConnection(origin, link2))
joint1 = EqualityConstraint(Spherical(origin, link1,zeros(3),zeros(3)))
joint2 = EqualityConstraint(Spherical(link1, link2, [0;0;2.],zeros(3)))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box1;box2]


mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)
setPosition!(mech,link1,q = Quaternion(RotZ(pi/2)))
setPosition!(mech,link2,x = [0.;0;2.],q = Quaternion(RotY(pi/2)))
setVelocity!(mech,link1,ω = [1.;0;0])
setVelocity!(mech,link1,link2, Δω = [10.;0;0.], p1=[0;0;1.], p2 = [0;0;-1.])


storage = simulate!(mech, 10., record = true)
visualize!(mech, storage, shapes)
