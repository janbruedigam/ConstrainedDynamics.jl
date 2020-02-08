using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width,depth = 0.5, 0.5
box1 = Box(width,depth,length1,length1,color=RGBA(1.,1.,0.))
box2 = Box(width,depth,length1,length1,color=RGBA(1.,0.,0.))

# Links
origin = Origin{Float64}()

link1 = Body(box1)
setInit!(origin,link1,[0.;1.;0.],zeros(3))

link2 = Body(box2)
setInit!(origin,link2,[0.;-1.;0.],zeros(3))

# Constraints
joint1 = Constraint(OriginConnection(origin,link1))
joint2 = Constraint(OriginConnection(origin,link2))

links = [link1;link2]
constraints = [joint1;joint2]
shapes = [box1;box2]


mech = Mechanism(origin, links, constraints,g=-5.)
link1.x[2] = [-0.05;0.95;0.05]
link2.x[2] = [-0.025;-0.975;0.05]
link2.q[2] = Quaternion(RotZ(0.1))

simulate!(mech,save=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
