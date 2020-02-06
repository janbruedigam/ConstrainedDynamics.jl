using Rotations

!(@isdefined FullCordDynamics) && include(joinpath("..", "FullCordDynamics.jl"))
using Main.FullCordDynamics

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width,depth = 0.1, 0.1
box = Box(width,depth,length1,length1)

p1 = [0.0;0.0;length1/2] # joint connection point

# Initial orientation
ϕ1 = π/2
q1 = Quaternion(RotX(ϕ1))

# Links
origin = Origin{Float64}()

link1 = Link(box)
setInit!(origin,link1,zeros(3),p1,q=q1)

# Constraints
joint_between_origin_and_link1 = Constraint(Socket(origin,link1,zeros(3),p1),Axis(origin,link1,joint_axis))

links = [link1]
constraints = [joint_between_origin_and_link1]
shapes = [box]


robot = Robot(origin, links, constraints)

simulate!(robot,save=true)
FullCordDynamics.visualize(robot,shapes)
