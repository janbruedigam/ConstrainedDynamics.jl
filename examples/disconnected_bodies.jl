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
joint1 = EqualityConstraint(OriginConnection(origin, link1))
joint2 = EqualityConstraint(OriginConnection(origin, link2))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box1;box2]


mech = Mechanism(origin, links, constraints, g = -5., shapes = shapes)
setPosition!(mech,link1,x = [0.;1.;0.])
setPosition!(mech,link2,x = [0.;-1.;0.])
setVelocity!(mech,link1,v = [-5.;-5;5])
setVelocity!(mech,link2,v = [-2.5;2.5;5],ω = [2;4;10.])


simulate!(mech,save = true)
visualize!(mech)
