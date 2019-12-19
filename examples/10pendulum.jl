using Rotations
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 0.5
x,y = .05,.05
b1 = Box(x,y,l1,l1,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11
vert1 = [[vert11];[vert12]]

# Initial orientation
phi = .7
q1=q2=q3=q4=q5=q6=q7=q8=q9=q10 = Quaternion(RotX(phi))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

link2 = Link(b1)
setInit!(link1,link2,vert12,vert11,q=q2)

link3 = Link(b1)
setInit!(link2,link3,vert12,vert11,q=q3)

link4 = Link(b1)
setInit!(link3,link4,vert12,vert11,q=q4)

link5 = Link(b1)
setInit!(link4,link5,vert12,vert11,q=q5)

link6 = Link(b1)
setInit!(link5,link6,vert12,vert11,q=q6)

link7 = Link(b1)
setInit!(link6,link7,vert12,vert11,q=q7)

link8 = Link(b1)
setInit!(link7,link8,vert12,vert11,q=q8)

link9 = Link(b1)
setInit!(link8,link9,vert12,vert11,q=q9)

link10 = Link(b1)
setInit!(link9,link10,vert12,vert11,q=q10)

# Constraints
jointb1 = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex))
joint12 = Constraint(Socket(link1,link2,vert12,vert11),Axis(link1,link2,ex))
joint23 = Constraint(Socket(link2,link3,vert12,vert11),Axis(link2,link3,ex))
joint34 = Constraint(Socket(link3,link4,vert12,vert11),Axis(link3,link4,ex))
joint45 = Constraint(Socket(link4,link5,vert12,vert11),Axis(link4,link5,ex))
joint56 = Constraint(Socket(link5,link6,vert12,vert11),Axis(link5,link6,ex))
joint67 = Constraint(Socket(link6,link7,vert12,vert11),Axis(link6,link7,ex))
joint78 = Constraint(Socket(link7,link8,vert12,vert11),Axis(link7,link8,ex))
joint89 = Constraint(Socket(link8,link9,vert12,vert11),Axis(link8,link9,ex))
joint910 = Constraint(Socket(link9,link10,vert12,vert11),Axis(link9,link10,ex))

links = [link1;link2;link3;link4;link5;link6;link7;link8;link9;link10]
constraints = [jointb1;joint12;joint23;joint34;joint45;joint56;joint67;joint78;joint89;joint910]
shapes = [b1]

bot = Robot(origin,links, constraints;tend=20.)

simulate!(bot,save=true,debug=false)
FullCordDynamics.visualize(bot,shapes)
