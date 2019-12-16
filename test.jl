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
origin = Link{Float64}()

link1 = Link{6}(m1,J1,vert1)
setInit!(origin,link1,[1;1],q=q1)

link2 = Link{6}(m2,J2,vert2)
setInit!(link1,link2,[2;1],q=q2)

link3 = Link{6}(m1,J1,vert1)
setInit!(origin,link3,[1;1],q=q3)

link4 = Link{6}(m2,J2,vert2)
setInit!(link3,link4,[2;1],q=q4)

# Constraints
joint0to1 = Combined2(Socket(origin,link1,1,1),Axis(origin,link1,ex))
joint1to23 = Combined3(Socket(link1,link2,2,1),Axis(link1,link2,ex),SocketYZ(link1,link3,1,1))
joint1to2 = Combined2(Socket(link1,link2,2,1),Axis(link1,link2,ex))
joint1to3 = Combined2(Socket(link1,link3,1,1),Axis(link1,link3,ex))

joint3to4 = Combined2(Socket(link3,link4,2,1),Axis(link3,link4,ex))
joint2to4 = Combined2(Socket(link2,link4,2,2),Axis(link2,link4,ex))


links = [link1; link2; link3; link4]
constraints = [joint0to1; joint1to23; joint3to4; joint2to4]
# constraints = [joint0to1; joint1to23; joint3to4]
# constraints = [joint0to1; joint1to2; joint1to3; joint3to4]
# links = [link1]
# constraints = [joint0to1]


bot = Robot(origin,links, constraints)

sim!(bot,save=true,debug=false)
include(joinpath("util", "visualize.jl"))
