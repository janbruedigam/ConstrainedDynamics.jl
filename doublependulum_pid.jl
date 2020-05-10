using ConstrainedDynamics
using ConstrainedControl

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1
box = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))

p2 = [0.0;0.0;length1 / 2] # joint connection point

# Initial orientation
ϕ1 = 0
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(box)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, zeros(3), p2, joint_axis))
joint2 = EqualityConstraint(Revolute(link1, link2, -p2, p2, joint_axis))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box]


mech = Mechanism(origin, links, constraints, shapes = shapes, tend = 10.)
setPosition!(mech,origin,link1,p2 = p2,Δq = q1)
setPosition!(mech,origin,link2,p1=-p2,p2 = p2,Δq = q1)

pid = PID(mech, getfield.(constraints,:id), [pi/2;-pi/4], P = [10.;0.], I = [10.;0], D = [5.;0.])


simulate!(mech,pid,save = true)
visualize!(mech)