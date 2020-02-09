using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
# include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width,depth = 0.5, 0.5
box1 = Box(width,depth,length1,length1,color=RGBA(1.,1.,0.))

# Links
origin = Origin{Float64}()

link1 = Body(box1)
setInit!(origin,link1,[0.;0.;2.],zeros(3))

# Constraints
joint1 = Constraint(OriginConnection(origin,link1))

links = [link1]
constraints = [joint1]
shapes = [box1]


mech = Mechanism(origin, links, constraints,g=-9.81,tend=10.)
link1.x[2]=[0.01;0.0;2.]

simulate_ip!(mech,save=true)

MaximalCoordinateDynamics.visualize(mech,shapes)
