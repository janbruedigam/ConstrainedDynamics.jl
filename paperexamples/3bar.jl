using Rotations
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
l2 = sqrt(2)/2
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))
b2 = Box(x,y,l2,l2,color=RGBA(1.,1.,0.))
b3 = Box(x,y,l1,l1,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

vert21 = [0.;0.;l2/2]
vert22 = -vert21

# Initial orientation
phi1, phi2, phi3 = 2.842889445268244, 2.842889445268244-1.209429202890086, pi
q1, q2, q3 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2)), Quaternion(RotX(phi3))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

link2 = Link(b2)
setInit!(link1,link2,vert12,vert21,q=q2)

link3 = Link(b3)
setInit!(origin,link3,[0;1.;0],vert11,q=q3)

# Constraints
joint0to13 = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex),SocketYZ(origin,link3,[0;1.;0],vert11))
joint1to2 = Constraint(Socket(link1,link2,vert12,vert21),Axis(link1,link2,ex))
joint2to3 = Constraint(Socket(link2,link3,vert22,vert12),Axis(link2,link3,ex))


links = [link1; link2; link3]
constraints = [joint0to13; joint1to2; joint2to3]
shapes = [b1,b2,b3]

bot = Robot(origin,links, constraints,tend=600.)

drift2 = simulate!(bot,save=true,debug=false)
FullCordDynamics.visualize(bot,shapes)
