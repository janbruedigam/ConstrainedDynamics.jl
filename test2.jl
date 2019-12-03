using Rotations
using BenchmarkTools
using TimerOutputs
using Plots

if !@isdefined includeFlag
    include("FullCordDynamics.jl")
    includeFlag = true
end
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
baseLink = Link(Float64,0,Np=2)

link1 = Link(m1,J1,vert1,1)
x1 =  initialPosition(baseLink,link1,q1,[1;1])
setInit!(link1,x=x1,q=q1)

link2 = Link(m1,J1,vert1,1)
x2 =  initialPosition(link1,link2,q2,[2;1])
setInit!(link2,x=x2,q=q2)

link3 = Link(m1,J1,vert1,1)
x3 =  initialPosition(link2,link3,q3,[2;1])
setInit!(link3,x=x3,q=q3)

link4 = Link(m1,J1,vert1,1)
x4 =  initialPosition(link3,link4,q4,[2;1])
setInit!(link4,x=x4,q=q4)

link5 = Link(m1,J1,vert1,1)
x5 =  initialPosition(link4,link5,q5,[2;1])
setInit!(link5,x=x5,q=q5)

link6 = Link(m1,J1,vert1,1)
x6 =  initialPosition(link5,link6,q6,[2;1])
setInit!(link6,x=x6,q=q6)

link7 = Link(m1,J1,vert1,1)
x7 =  initialPosition(link6,link7,q7,[2;1])
setInit!(link7,x=x7,q=q7)

link8 = Link(m1,J1,vert1,1)
x8 =  initialPosition(link7,link8,q8,[2;1])
setInit!(link8,x=x8,q=q8)

link9 = Link(m1,J1,vert1,1)
x9 =  initialPosition(link8,link9,q9,[2;1])
setInit!(link9,x=x9,q=q9)

link10 = Link(m1,J1,vert1,1)
x10 =  initialPosition(link9,link10,q10,[2;1])
setInit!(link10,x=x10,q=q10)

# Constraints
jointb = Combined(FixedPosition(baseLink,1),FixedOrientation(baseLink))
jointb1 = Combined(Socket(baseLink,link1,1,1),Axis(baseLink,link1,ex))
joint12 = Combined(Socket(link1,link2,2,1),Axis(link1,link2,ex))
joint23 = Combined(Socket(link2,link3,2,1),Axis(link2,link3,ex))
joint34 = Combined(Socket(link3,link4,2,1),Axis(link3,link4,ex))
joint45 = Combined(Socket(link4,link5,2,1),Axis(link4,link5,ex))
joint56 = Combined(Socket(link5,link6,2,1),Axis(link5,link6,ex))
joint67 = Combined(Socket(link6,link7,2,1),Axis(link6,link7,ex))
joint78 = Combined(Socket(link7,link8,2,1),Axis(link7,link8,ex))
joint89 = Combined(Socket(link8,link9,2,1),Axis(link8,link9,ex))
joint910 = Combined(Socket(link9,link10,2,1),Axis(link9,link10,ex))

links = [baseLink; link1;link2;link3;link4;link5;link6;link7;link8;link9;link10]
constraints = [jointb; jointb1;joint12;joint23;joint34;joint45;joint56;joint67;joint78;joint89;joint910]

# links = [baseLink; link1]
# constraints = [jointb; jointb1]
# links = [baseLink]
# constraints = [jointb]


bot = Robot(links, constraints, root=length(links)+1)

# sim!(bot)
