using ConstrainedDynamics
using ConstrainedDynamicsVis

joint_axis = [1.0;0.0;0.0]
l = 1.0
w, d = .1, .1
b1 = Box(w, d, l, l)
b2 = Box(w, d, l, l)

vert11 = [0.;0.;l / 2]
vert12 = -vert11
vert21 = [0.;0.;l / 2]

origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b2)
joint0to1 = EqualityConstraint(Fixed(origin, link1; p2=vert11))
joint1to2 = EqualityConstraint(Revolute(link1, link2,joint_axis; p1=vert12, p2=vert21))

links = [link1;link2]
constraints = [joint0to1;joint1to2]
shapes = [b1,b2]

mech = Mechanism(origin, links, constraints, shapes=shapes)
setPosition!(origin,link1,p2 = [1;0;0])
setPosition!(link1,link2,p1 = [1;2;3],p2 = [1;2;4],Δq = UnitQuaternion(RotX(pi/4)))

initializeConstraints!(mech,[1]) 

storage = simulate!(mech, 10, record = true)
visualize(mech, storage, shapes)