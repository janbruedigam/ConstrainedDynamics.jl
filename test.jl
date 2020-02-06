using Rotations
using BenchmarkTools
using Plots

(@isdefined MaximalCoordinateDynamics) ? nothing : include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
l2 = sqrt(2)/2
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))
b2 = Box(x,y,l2,l2,color=RGBA(1.,1.,0.))
b3 = Box(x,y,l1,l1,color=RGBA(1.,0.,0.))
b4 = Box(x,y,l2,l2,color=RGBA(1.,0.,0.))
b5 = Box(0.1,0.1,0.1,0.1,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

vert21 = [0.;0.;l2/2]
vert22 = -vert21

# Initial orientation
phi1, phi2, phi3, phi4 = 0.4,0.,0.,0.#pi/2, -pi/4, 0., 3*pi/4.
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
setInit!(origin,link5,zeros(3),zeros(3))
setInit!(link5,link1,zeros(3),vert11,q=q1)
setInit!(link1,link2,vert12,vert21,q=q2)

# Constraints

# joint0to1 = Constraint(Socket(origin,link1,zeros(3),vert11))
# joint1to23 = Constraint(Socket(link1,link2,vert12,vert21),SocketYZ(link1,link3,vert11,vert11))
#
# joint3to4 = Constraint(Socket(link3,link4,vert12,vert21))
# joint2to4 = Constraint(Socket(link2,link4,vert22,vert22))
#
# joint23to4 = Constraint(Socket(link3,link4,vert12,vert21))
# joint333to4 = Constraint(Socket(link2,link4,vert22,vert22))
#
# links = [link1; link2; link3; link4]
# constraints = [joint0to1; joint1to23; joint3to4; joint2to4]
shapes = [b1,b2,b3,b4,b5]

oc = Constraint(OriginConnection(origin,link1))
# testjoint = Constraint(MatchedOrientation(link1,link2,q1))
# testjoint = Constraint(MatchedOrientation(origin,link1))
# testjoint = Constraint(Line(origin,link1,[0.;0.;1.]))
testjoint = Constraint(MatchedOrientation(origin,link5),Line(origin,link5,[0.;1.;0.],zeros(3),zeros(3)))
joint5to1 = Constraint(Socket(link5,link1,zeros(3),vert11))
joint1to2 = Constraint(Socket(link1,link2,vert12,vert21))

# bot = Robot(origin,links,constraints)
# bot = Robot(origin,[link1;link2],[oc;testjoint])
bot = Robot(origin,[link5;link1;link2],[testjoint;joint5to1;joint1to2])

simulate!(bot,save=true,debug=false)
MaximalCoordinateDynamics.visualize(bot,shapes)
