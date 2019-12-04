using Rotations
using BenchmarkTools
using TimerOutputs
using Plots

(@isdefined FullCordDynamics) ? nothing : include("FullCordDynamics.jl")
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.
m1, J1 = box(.1,.1,l1,l1)
vert11 = [0.;0.;l1/2]
vert12 = -vert11
vert1 = [[vert11];[vert12]]

# Initial orientation
phi = .2
q1=q2=q3=q4=q5=q6=q7=q8=q9=q10 = Quaternion(RotX(phi))

# Links
origin = Link(Float64)

link1 = Link(m1,J1,vert1,1)
setInit!(link1,origin,[1;1],q=q1)

link2 = Link(m1,J1,vert1,1)
setInit!(link2,link1,[2;1],q=q2)

link3 = Link(m1,J1,vert1,1)
setInit!(link3,link2,[2;1],q=q3)

link4 = Link(m1,J1,vert1,1)
setInit!(link4,link3,[2;1],q=q4)

link5 = Link(m1,J1,vert1,1)
setInit!(link5,link4,[2;1],q=q5)

link6 = Link(m1,J1,vert1,1)
setInit!(link6,link5,[2;1],q=q6)

link7 = Link(m1,J1,vert1,1)
setInit!(link7,link6,[2;1],q=q7)

link8 = Link(m1,J1,vert1,1)
setInit!(link8,link7,[2;1],q=q8)

link9 = Link(m1,J1,vert1,1)
setInit!(link9,link8,[2;1],q=q9)

link10 = Link(m1,J1,vert1,1)
setInit!(link10,link9,[2;1],q=q10)

# Constraints
jointb1 = Combined(Socket(origin,link1,1,1),Axis(origin,link1,ex))
joint12 = Combined(Socket(link1,link2,2,1),Axis(link1,link2,ex))
joint23 = Combined(Socket(link2,link3,2,1),Axis(link2,link3,ex))
joint34 = Combined(Socket(link3,link4,2,1),Axis(link3,link4,ex))
joint45 = Combined(Socket(link4,link5,2,1),Axis(link4,link5,ex))
joint56 = Combined(Socket(link5,link6,2,1),Axis(link5,link6,ex))
joint67 = Combined(Socket(link6,link7,2,1),Axis(link6,link7,ex))
joint78 = Combined(Socket(link7,link8,2,1),Axis(link7,link8,ex))
joint89 = Combined(Socket(link8,link9,2,1),Axis(link8,link9,ex))
joint910 = Combined(Socket(link9,link10,2,1),Axis(link9,link10,ex))

jointb2 = Combined(Socket(origin,link2,1,1),Axis(origin,link2,ex))
jointb3 = Combined(Socket(origin,link3,1,1),Axis(origin,link3,ex))
jointb4 = Combined(Socket(origin,link4,1,1),Axis(origin,link4,ex))
jointb5 = Combined(Socket(origin,link5,1,1),Axis(origin,link5,ex))
jointb6 = Combined(Socket(origin,link6,1,1),Axis(origin,link6,ex))
jointb7 = Combined(Socket(origin,link7,1,1),Axis(origin,link7,ex))
jointb8 = Combined(Socket(origin,link8,1,1),Axis(origin,link8,ex))
jointb9 = Combined(Socket(origin,link9,1,1),Axis(origin,link9,ex))
jointb10 = Combined(Socket(origin,link10,1,1),Axis(origin,link10,ex))

links = [link1;link2;link3;link4;link5;link6;link7;link8;link9;link10]
constraints = [jointb1;joint12;joint23;joint34;joint45;joint56;joint67;joint78;joint89;joint910]
# constraints = [jointb1;jointb2;jointb3;jointb4;jointb5;jointb6;jointb7;jointb8;jointb9;jointb10]

# links = [baseLink; link1;link2;link3;link4;link5]
# constraints = [jointb1;joint12;joint23;joint34;joint45]

# links = [baseLink; link1]
# constraints = [jointb1]


bot = Robot(origin,links, constraints)

# sim!(bot,save=true)
# trajS = trajSFunc(bot)
