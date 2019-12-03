using Rotations
using BenchmarkTools
using TimerOutputs
using Plots

if !@isdefined includeFlag
    include("FullCordDynamics.jl")
    includeFlag = true
end
using Main.FullCordDynamics

# const to = TimerOutput()

# Parameters
ex = [1.;0.;0.]

l1 = 1.
m1, J1 = box(.1,.1,l1,l1)
vert11 = [0.;0.;l1/2]
vert12 = -vert11
vert1 = [[vert11];[vert12]]

l2 = sqrt(2)/2
m2, J2 = box(.1,.1,l2,l2)
vert21 = [0.;0.;l2/2]
vert22 = -vert21
vert2 = [[vert21];[vert22]]

# Initial orientation
phi1, phi2, phi3, phi4 = pi/2, -pi/4, 0., 3*pi/4.
q1, q2, q3, q4 = Quaternion(RotX(phi1)), Quaternion(RotX(phi2)), Quaternion(RotX(phi3)), Quaternion(RotX(phi4))

# Links
baseLink = Link(Float64,0,Np=2)

link1 = Link(m1,J1,vert1,1)
x1 =  initialPosition(baseLink,link1,q1,[1;1])
setInit!(link1,x=x1,q=q1)

link2 = Link(m2,J2,vert2,1)
x2 = initialPosition(link1,link2,q2,[2;1])
setInit!(link2,x=x2,q=q2)

link3 = Link(m1,J1,vert1,1)
x3 = initialPosition(baseLink,link3,q3,[1;1])
setInit!(link3,x=x3,q=q3)

link4 = Link(m2,J2,vert2,1)
x4 = initialPosition(link3,link4,q4,[2;1])
setInit!(link4,x=x4,q=q4)

# Constraints
jointb = Combined(FixedPosition(baseLink,1),FixedOrientation(baseLink))
jointb1 = Combined(Socket(baseLink,link1,1,1),Axis(baseLink,link1,ex))
joint12 = Combined(Socket(link1,link2,2,1),Axis(link1,link2,ex))
jointb3 = Combined(Socket(baseLink,link3,1,1),Axis(baseLink,link3,ex))
joint34 = Combined(Socket(link3,link4,2,1),Axis(link3,link4,ex))

links = [baseLink; link1; link2; link3; link4]
constraints = [jointb; jointb1; joint12; jointb3; joint34]
# links = [baseLink; link1]
# constraints = [jointb; jointb1]


bot = Robot(links, constraints, root=length(links)+1)

sim!(bot,save=true)
trajS = trajSFunc(bot)
include("visualizeTwoTwoBar.jl")
