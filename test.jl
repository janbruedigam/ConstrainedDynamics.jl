using Rotations
using BenchmarkTools
using TimerOutputs
using Plots

(@isdefined FullCordDynamics) ? nothing : include("FullCordDynamics.jl")
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
origin = Link(Float64)

link1 = Link(m1,J1,vert1,1)
setInit!(link1,origin,[1;1],q=q1)

link2 = Link(m2,J2,vert2,1)
setInit!(link2,link1,[2;1],q=q2)

link3 = Link(m1,J1,vert1,1)
setInit!(link3,origin,[1;1],q=q3)

link4 = Link(m2,J2,vert2,1)
setInit!(link4,link3,[2;1],q=q4)

# Constraints
jointb1 = Combined(Socket(origin,link1,1,1),Axis(origin,link1,ex))
joint12 = Combined(Socket(link1,link2,2,1),Axis(link1,link2,ex))
jointb3 = Combined(Socket(origin,link3,1,1),Axis(origin,link3,ex))
joint34 = Combined(Socket(link3,link4,2,1),Axis(link3,link4,ex))

# links = [link1; link2; link3; link4]
# constraints = [jointb1; joint12; jointb3; joint34]
links = [link1]
constraints = [jointb1]


bot = Robot(origin,links, constraints)

# sim!(bot,save=true)
# trajS = trajSFunc(bot)
# include("visualizeTwoTwoBar.jl")
