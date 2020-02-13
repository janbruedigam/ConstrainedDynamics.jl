using Rotations
using Plots: RGBA
using StaticArrays

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

link1 = Body(box1)
setInit!(origin,link1,[0.;-2;1.5],zeros(3))

# link1 = Body(b1)
# setInit!(origin,link1,[0;0;0.0],[0;0;-length1],q = Quaternion(RotX(0.3)))
#
# link2 = Body(b1)
# setInit!(link1,link2,[0;0;length1],[0;0;length1],q = Quaternion(RotX(-0.3)))
#
link2 = Body(b2)
setInit!(link1,link2,[length1/2;length1/2;-length1/2],zeros(3))

link3 = Body(b2)
setInit!(link1,link3,[length1/2;-length1/2;-length1/2],zeros(3))

link4 = Body(b2)
setInit!(link1,link4,[-length1/2;length1/2;-length1/2],zeros(3))

link5 = Body(b2)
setInit!(link1,link5,[-length1/2;-length1/2;-length1/2],zeros(3))

link6 = Body(b2)
setInit!(link1,link6,[length1/2;length1/2;length1/2],zeros(3))

link7 = Body(b2)
setInit!(link1,link7,[length1/2;-length1/2;length1/2],zeros(3))

link8 = Body(b2)
setInit!(link1,link8,[-length1/2;length1/2;length1/2],zeros(3))

link9 = Body(b2)
setInit!(link1,link9,[-length1/2;-length1/2;length1/2],zeros(3))

# # Constraints
joint1 = InequalityConstraint(link1)
joint2 = InequalityConstraint(link2)
joint3 = InequalityConstraint(link3)
joint4 = InequalityConstraint(link4)
joint5 = InequalityConstraint(link5)
joint6 = InequalityConstraint(link6)
joint7 = InequalityConstraint(link7)
joint8 = InequalityConstraint(link8)
joint9 = InequalityConstraint(link9)
#
# links = [link1]
# constraints = [joint1]
# ineqs = [joint2]
# # ineqs = InequalityConstraint{Float64}[]
# shapes = [box1]

# Constraints
# joint0to1 = EqualityConstraint(Revolute(origin,link1,[0.;0;2.],zeros(3),joint_axis))
# joint0to1 = EqualityConstraint(OriginConnection(origin,link1))
# joint1to2 = EqualityConstraint(Revolute(link1,link2,[0.;0;length1],[0.;0;length1],joint_axis))
# joint1to3 = EqualityConstraint(Fixed(link1,link3,[0.;0;-length1],zeros(3)))
# joint2to4 = EqualityConstraint(Fixed(link2,link4,[0.;0;-length1],zeros(3)))
# joint3ineq = InequalityConstraint(link3)
# joint4ineq = InequalityConstraint(link4)

# links = [link1;link2;link3;link4]
# constraints = [joint0to1;joint1to2;joint1to3;joint2to4]
# ineqs = [joint3ineq;joint4ineq]
# ineqs = InequalityConstraint{Float64}[]

# joint0to1 = EqualityConstraint(Revolute(origin,link1,[0.;0;1.5],[0;2;0.],joint_axis))
joint0to1 = EqualityConstraint(OriginConnection(origin,link1))
joint1to2 = EqualityConstraint(Fixed(link1,link2,[length1/2;length1/2;-length1/2],zeros(3)))
joint1to3 = EqualityConstraint(Fixed(link1,link3,[length1/2;-length1/2;-length1/2],zeros(3)))
joint1to4 = EqualityConstraint(Fixed(link1,link4,[-length1/2;length1/2;-length1/2],zeros(3)))
joint1to5 = EqualityConstraint(Fixed(link1,link5,[-length1/2;-length1/2;-length1/2],zeros(3)))
joint1to6 = EqualityConstraint(Fixed(link1,link6,[length1/2;length1/2;length1/2],zeros(3)))
joint1to7 = EqualityConstraint(Fixed(link1,link7,[length1/2;-length1/2;length1/2],zeros(3)))
joint1to8 = EqualityConstraint(Fixed(link1,link8,[-length1/2;length1/2;length1/2],zeros(3)))
joint1to9 = EqualityConstraint(Fixed(link1,link9,[-length1/2;-length1/2;length1/2],zeros(3)))
links = [link1;link2;link3;link4;link5;link6;link7;link8;link9]
constraints = [joint0to1;joint1to2;joint1to3;joint1to4;joint1to5;joint1to6;joint1to7;joint1to8;joint1to9]
ineqs = [joint2;joint3;joint4;joint5;joint6;joint7;joint8;joint9]
shapes = [box1;b1;b2]


mech = Mechanism(origin, links,constraints,ineqs,g=-9.81,tend=10.)
link1.q[2] = Quaternion(AngleAxis(-rand()-0.2,rand(3)...))
# link1.q[2] = Quaternion(SVector([0.885818;-0.0789202;-0.274472;-0.365735]...))
qtemp = link1.q[2]
# 0.9127362490430289
# -0.19667611676218716
# -0.01932481073253633
# -0.3575718060311244
link1.x[2] = [0.;-1.95;1.55]

simulate_ip!(mech,save=true,debug=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
