using Rotations
using Plots

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))

vert11 = [0.;0.;l1/2]

# Initial orientation
phi1 = pi/2
q1 = Quaternion(RotX(phi1))

# Links
origin = Origin{Float64}()

link1 = Link(b1)
setInit!(origin,link1,zeros(3),vert11,q=q1)

# Constraints
joint0to1 = Constraint(Socket(origin,link1,zeros(3),vert11),Axis(origin,link1,ex))

links = [link1]
constraints = [joint0to1]
shapes = [b1]


bot = Robot(origin,links, constraints)

simulate!(bot,save=true,debug=false)
FullCordDynamics.visualize(bot,shapes)
