using Rotations
using BenchmarkTools
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))
b2 = Box(0.1,0.1,0.1,1.,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

# Initial orientation
phi1, phi2 = pi/4, -pi/4
q1, q2 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2))

# Links
origin = Origin{Float64}()

shifter = Link(b2)
setInit!(origin,shifter,[0.;-.5*sqrt(2);0.],zeros(3))

link1 = Link(b1)
setInit!(shifter,link1,zeros(3),vert11,q=q1)

link2 = Link(b1)
setInit!(origin,link2,zeros(3),vert11,q=q2)

link3 = Link(b1)
setInit!(link1,link3,vert12,vert11,q=q2)

link4 = Link(b1)
setInit!(link2,link4,vert12,vert11,q=q1)

link5 = Link(b1)
setInit!(link3,link5,vert12,vert11,q=q1)

link6 = Link(b1)
setInit!(link4,link6,vert12,vert11,q=q2)

# Constraints

shapes = [b1,b2]

joint0tos2 = Constraint(MatchedOrientation(origin,shifter),Line(origin,shifter,[0.;1.;0.],zeros(3),zeros(3)),SocketYZ(origin,link2,zeros(3),vert11))
jointsto1 = Constraint(SocketYZ(shifter,link1,zeros(3),vert11))
joint1to23 = Constraint(SocketYZ(link1,link2,zeros(3),zeros(3)))#,SocketYZ(link1,link3,vert12,vert11))
joint2to4 = Constraint(SocketYZ(link2,link4,vert12,vert11))
joint3to45 = Constraint(SocketYZ(link3,link4,zeros(3),zeros(3)),SocketYZ(link3,link5,vert12,vert11))
joint4to6 = Constraint(SocketYZ(link4,link6,vert12,vert11))
joint5to6 = Constraint(SocketYZ(link5,link6,zeros(3),zeros(3)))

links = [shifter;link1;link2]
constraints = [joint0tos2;jointsto1;joint1to23]

bot = Robot(origin,links,constraints)

simulate!(bot,save=true,debug=false)
FullCordDynamics.visualize(bot,shapes)
