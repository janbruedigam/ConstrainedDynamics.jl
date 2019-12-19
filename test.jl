using Rotations
using BenchmarkTools
using TimerOutputs
using Plots

(@isdefined FullCordDynamics) ? nothing : include("FullCordDynamics.jl")
using Main.FullCordDynamics

# const to = TimerOutput()

# Parameters
ex = [1.;0.;0.]

l1 = 1.0*0.3
l2 = sqrt(2)/2
x,y = .1,.1
# b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))
b1 = Box(l1,y,y,l1,color=RGBA(1.,1.,0.))
b2 = Box(x,y,l2,l2,color=RGBA(1.,1.,0.))
b3 = Box(x,y,l1,l1,color=RGBA(1.,0.,0.))
b4 = Box(x,y,l2,l2,color=RGBA(1.,0.,0.))
b5 = Box(.1,1.,1.,2.,color=RGBA(1.,0.,0.))

# vert11 = [0.;0.;l1/2]
# vert12 = -vert11
vert11 = [l1/2;0.;0.]
vert12 = -vert11

vert21 = [0.;0.;l2/2]
vert22 = -vert21

# Initial orientation
phi1, phi2, phi3, phi4 = 0*pi/2, -pi/4, 0., 3*pi/4.
q1, q2, q3, q4 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2)), Quaternion(RotX(phi3)), Quaternion(RotX(phi4))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

link2 = Link(b2)
setInit!(link1,link2,vert12,vert21,q=q2)

link3 = Link(b3)
setInit!(link1,link3,vert11,vert11,q=q3)

link4 = Link(b4)
setInit!(link3,link4,vert12,vert21,q=q4)

link5 = Link(b5)
setInit!(link1,link5,vert12,zeros(3),q=q1,τ=[0.0;0.;0.])

# Constraints
socket0to1 = Constraint(Socket(origin,link1,zeros(3),vert11))
socket1to23 = Constraint(Socket(link1,link2,vert12,vert21),SocketYZ(link1,link3,vert11,vert11))
socket3to4 = Constraint(Socket(link3,link4,vert12,vert21))
socket2to4 = Constraint(Socket(link2,link4,vert22,vert22))

joint0to1 = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex))
joint1to23 = Constraint(Socket(link1,link2,vert12,vert21),Axis(link1,link2,ex),SocketYZ(link1,link3,vert11,vert11))
joint1to2 = Constraint(Socket(link1,link2,vert12,vert21),Axis(link1,link2,ex))
joint1to3 = Constraint(Socket(link1,link3,vert11,vert11),Axis(link1,link3,ex))

joint3to4 = Constraint(Socket(link3,link4,vert12,vert21),Axis(link3,link4,ex))
joint2to4 = Constraint(Socket(link2,link4,vert22,vert22),Axis(link2,link4,ex))

joint1to5= Constraint(Socket(link1,link5,vert12,zeros(3)),Axis(link1,link5,ex))


# links = [link1; link2; link3; link4]
# constraints = [joint0to1; joint1to23; joint3to4; joint2to4]
# constraints = [joint0to1; joint1to23; joint3to4]
# constraints = [joint0to1; joint1to2; joint1to3; joint3to4]
# constraints = [socket0to1; joint1to23; joint3to4; joint2to4]
# constraints = [socket0to1;socket1to23;socket3to4;socket2to4]
links = [link1;link5]
constraints = [socket0to1;joint1to5]
shapes = [b1,b2,b3,b4,b5]


bot = Robot(origin,links, constraints;g=-9.81,dt=0.001)
link5.q[2] = Quaternion(RotX(0.05))

sim!(bot,save=true,debug=false)
FullCordDynamics.visualize(bot,shapes)
