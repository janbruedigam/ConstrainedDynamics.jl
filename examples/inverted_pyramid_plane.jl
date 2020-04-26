using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width, depth = 0.5, 0.5
box1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))

# Links
origin = Origin{Float64}()
link1 = Body(box1)

# Constraints
joint1 = InequalityConstraint(Impact(link1, [0;-0.1;1.0]), Impact(link1, [0;0.1;1.0]))
joint2 = InequalityConstraint(Impact(link1, [0.1;0;1.0]))
joint3 = InequalityConstraint(Impact(link1, [-0.1;0;1.0]))

links = [link1]
ineqs = [joint1;joint2;joint3]
shapes = [box1]


mech = Mechanism(origin, links, ineqs, shapes = shapes)
setVelocity!(mech,link1,v = [1;0.5;5])


simulate!(mech,save = true)
visualize!(mech)
