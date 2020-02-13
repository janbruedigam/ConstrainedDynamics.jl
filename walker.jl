using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
# include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width,depth = 0.5, 0.5
box1 = Box(width,depth,length1,1.,color=RGBA(1.,1.,0.))
b1 = Box(0.1,0.1,2length1,1.,color=RGBA(1.,1.,0.))
b2 = Box(0.1,0.1,.1,.1,color=RGBA(1.,0.,0.))

# Links
origin = Origin{Float64}()

# link1 = Body(box1)
# setInit!(origin,link1,[0.;0;100],zeros(3))

link1 = Body(b1)
setInit!(origin,link1,[0;0;1.],[0;0;-length1],q = Quaternion(RotX(0.3)))

link2 = Body(b1)
setInit!(link1,link2,[0;0;length1],[0;0;length1],q = Quaternion(RotX(-1.)))

link3 = Body(b2)
setInit!(link1,link3,[0.;0;-length1],zeros(3))

link4 = Body(b2)
setInit!(link2,link4,[0.;0;-length1],zeros(3))

# # Constraints
# joint1 = EqualityConstraint(OriginConnection(origin,link1))
# joint2 = InequalityConstraint(link1)
#
# links = [link1]
# constraints = [joint1]
# ineqs = [joint2]
# # ineqs = InequalityConstraint{Float64}[]
# shapes = [box1]

# Constraints
# joint0to1 = EqualityConstraint(Revolute(origin,link1,[0.;0;2.],zeros(3),joint_axis))
joint0to1 = EqualityConstraint(OriginConnection(origin,link1))
joint1to2 = EqualityConstraint(Revolute(link1,link2,[0.;0;length1],[0.;0;length1],joint_axis))
joint1to3 = EqualityConstraint(Fixed(link1,link3,[0.;0;-length1],zeros(3)))
joint2to4 = EqualityConstraint(Fixed(link2,link4,[0.;0;-length1],zeros(3)))
joint3ineq = InequalityConstraint(link3)
joint4ineq = InequalityConstraint(link4)

links = [link1;link2;link3;link4]
constraints = [joint0to1;joint1to2;joint1to3;joint2to4]
ineqs = [joint3ineq;joint4ineq]
# ineqs = InequalityConstraint{Float64}[]
# links = [link1]
# constraints = [joint1]
# ineqs = [joint2]
shapes = [box1;b1;b2]


mech = Mechanism(origin, links,constraints,ineqs,g=-9.81,tend=10.)
link1.q[2] = Quaternion(AngleAxis(-rand()*.1,rand(3)-ones(3)*0.5...))

simulate_ip!(mech,save=true,debug=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
