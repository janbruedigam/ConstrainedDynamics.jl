using ConstrainedDynamics


# Parameters
joint_axis = [.0;1.0;1.0]
v1 = -[0.5;0;0.]
v2 = -[0.;0.5;0.]
v3 = -[0.;0;0.5]

length1 = 1.
width, depth = 0.1, 0.1
box1 = Box(length1, width, depth, length1, color = RGBA(1., 1., 0.))
box2 = Box(width, length1, depth, length1, color = RGBA(1., 0., 0.))
box3 = Box(width, depth, length1, length1, color = RGBA(1., 0., 1.))

# Links
origin = Origin{Float64}()
link1 = Body(box1)
link2 = Body(box2)
link3 = Body(box3)

# Constraints
# joint1 = EqualityConstraint(Spherical(origin, link1,zeros(3),v1))
# joint2 = EqualityConstraint(Spherical(origin, link2,zeros(3),v2))
# joint3 = EqualityConstraint(Spherical(origin, link3,zeros(3),v3))
joint1 = EqualityConstraint(Revolute(origin, link1,zeros(3),v1,joint_axis))
joint2 = EqualityConstraint(Revolute(origin, link2,zeros(3),v2,joint_axis))
joint3 = EqualityConstraint(Revolute(origin, link3,zeros(3),v3,joint_axis))

links = [link1;link2;link3]
constraints = [joint1;joint2;joint3]
shapes = [box1;box2;box3]


mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)
setPosition!(origin,link1,Δq = Quaternion(RotZ(pi/2)),p2=v1)
setPosition!(origin,link2,Δq = Quaternion(RotX(0.)),p2=v2)
setPosition!(origin,link3,Δq = Quaternion(RotX(0.)),p2=v3)
# setVelocity!(origin,link1,Δω = [0.;1;1],p2=v1)
# setVelocity!(origin,link2,Δω = [0.;1;1],p2=v2)
# setVelocity!(origin,link3,Δω = [0.;1;1],p2=v3)
# setVelocity!(link2,v = [-2.5;2.5;5],ω = [2;4;10.])
setForce!(link1,τ=[0.1;0;0])

