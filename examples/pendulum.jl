using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1)

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = π / 2
q1 = UnitQuaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)

# Constraints
joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


mech = LinearMechanism(Mechanism(origin, links, constraints, shapes = shapes))
setPosition!(origin,link1,p2 = p2,Δq = q1)

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)
