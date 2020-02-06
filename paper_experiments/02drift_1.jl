using Rotations


!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
joint_axis = [1.;0.;0.]
m = 1. # mass
l1 = 1.0 # length1
l2 = sqrt(2)/2 # length2
x,y = .1,.1 # size of link
box1 = Box(x,y,l1,l1)
box2 = Box(x,y,l2,l2)
box3 = Box(x,y,l1,l1)

 # joint connection points
p1 = [0.;0.;l1/2]
p2 = -p1
p3 = [0.;0.;l2/2]
p4 = -p3

# Initial orientation
ϕ1, ϕ2, ϕ3 = 2.842889445268244, 2.842889445268244-1.209429202890086, pi
q1, q2, q3 = Quaternion(RotX(ϕ1)), Quaternion(RotX(ϕ2)), Quaternion(RotX(ϕ3))

# Links
origin = Origin{Float64}()

link1 = Link(box1)
setInit!(origin,link1,zeros(3),p1,q=q1)

link2 = Link(box2)
setInit!(link1,link2,p2,p3,q=q2)

link3 = Link(box3)
setInit!(origin,link3,[0;1.;0],p1,q=q3)

# Constraints
joint0to13 = Constraint(Socket(origin,link1,zeros(3),p1),Axis(origin,link1,joint_axis),SocketYZ(origin,link3,[0;1.;0],p1))
joint1to2 = Constraint(Socket(link1,link2,p2,p3),Axis(link1,link2,joint_axis))
joint2to3 = Constraint(Socket(link2,link3,p4,p2),Axis(link2,link3,joint_axis))


links = [link1; link2; link3]
constraints = [joint0to13; joint1to2; joint2to3]

# Mechanism
bot = Robot(origin,links, constraints,tend=600.)

drift = simulate_drift!(bot)
