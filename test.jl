using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
# include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.
width,depth = 0.1, 0.1
box1 = Box(width,depth,length1,length1,color=RGBA(1.,1.,0.))

# Links
origin = Origin{Float64}()

link1 = Body(box1)
setInit!(origin,link1,[0.;0;2],zeros(3))

link2 = Body(box1)
setInit!(origin,link2,[0.;length1/2;1.5],zeros(3),q=Quaternion(RotX(pi/2)))

# Constraints
joint1 = EqualityConstraint(Revolute(origin,link1,[0;0;2.5],[0;0;0.5],joint_axis))
joint2 = EqualityConstraint(Revolute(link1,link2,[0;0;-length1/2],[0;0;length1/2],joint_axis))
joint3 = InequalityConstraint(link2)

jointb = EqualityConstraint(OriginConnection(origin,link2))

links = [link1;link2]
constraints = [joint1;jointb]
ineqs = [joint3]
shapes = [box1]


mech = Mechanism(origin, links,constraints,ineqs,g=-9.81,tend=10.)

simulate_ip!(mech,save=true,debug=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
