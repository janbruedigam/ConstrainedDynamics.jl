using ConstrainedDynamics


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = π / 2
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, zeros(3), p2, joint_axis))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)

simulate!(mech,save = true)
visualize!(mech)
