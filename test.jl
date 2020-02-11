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
setInit!(origin,link1,[0.;0;2],zeros(3))

# Constraints
joint1 = EqualityConstraint(OriginConnection(origin,link1))
joint2 = InequalityConstraint(link1)

links = [link1]
constraints = [joint1]
ineqs = [joint2]
# ineqs = InequalityConstraint{Float64}[]
shapes = [box1]


mech = Mechanism(origin, links,constraints,ineqs,g=-9.81,tend=10.)
link1.x[2] = [0;0.01;2]

simulate_ip!(mech,save=true,debug=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
