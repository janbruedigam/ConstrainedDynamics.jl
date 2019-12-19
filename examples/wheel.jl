using Rotations
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0*0.3
d = .05
b1 = Box(l1,d,d,l1,color=RGBA(1.,1.,0.))
b5 = Box(.1,1.,1.,2.,color=RGBA(1.,0.,0.))

vert11 = [l1/2;0.;0.]
vert12 = -vert11

# Initial orientation
phi1 = 0.
q1 = Quaternion(RotX(phi1))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

link5 = Link(b5)
setInit!(link1,link5,vert12,zeros(3),q=q1,τ=[0.0;0.;0.])

# Constraints
socket0to1 = Constraint(Socket(origin,link1,zeros(3),vert11))
joint1to5= Constraint(Socket(link1,link5,vert12,zeros(3)),Axis(link1,link5,ex))

links = [link1;link5]
constraints = [socket0to1;joint1to5]
shapes = [b1;b5]


bot = Robot(origin,links, constraints;dt=0.001)
link5.q[2] = Quaternion(RotX(0.05))

simulate!(bot,save=true,debug=false)
FullCordDynamics.visualize(bot,shapes)
