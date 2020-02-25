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
setInit!(origin,link1,[0.;-2;1.5],zeros(3),F=[0;0.0;0])

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
cfr = .15
# joint1 = InequalityConstraint(Impact(link1,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint2 = InequalityConstraint(Impact(link2,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint3 = InequalityConstraint(Impact(link3,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint4 = InequalityConstraint(Impact(link4,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint5 = InequalityConstraint(Impact(link5,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint6 = InequalityConstraint(Impact(link6,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint7 = InequalityConstraint(Impact(link7,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint8 = InequalityConstraint(Impact(link8,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))
# joint9 = InequalityConstraint(Impact(link9,[0;-.1;01.0]),Impact(link1,[0;.1;01.0]))

joint1 = InequalityConstraint(Impact(link1,[0.0;.20;01.0]),Impact(link1,[0.0;-.20;01.0]))
joint2 = InequalityConstraint(Impact(link2,[0.0;.20;01.0]),Impact(link2,[0.0;-.20;01.0]))
joint3 = InequalityConstraint(Impact(link3,[0.0;.20;01.0]),Impact(link3,[0.0;-.20;01.0]))
joint4 = InequalityConstraint(Impact(link4,[0.0;.20;01.0]),Impact(link4,[0.0;-.20;01.0]))
joint5 = InequalityConstraint(Impact(link5,[0.0;.20;01.0]),Impact(link5,[0.0;-.20;01.0]))
joint6 = InequalityConstraint(Impact(link6,[0.0;.20;01.0]),Impact(link6,[0.0;-.20;01.0]))
joint7 = InequalityConstraint(Impact(link7,[0.0;.20;01.0]),Impact(link7,[0.0;-.20;01.0]))
joint8 = InequalityConstraint(Impact(link8,[0.0;.20;01.0]),Impact(link8,[0.0;-.20;01.0]))
joint9 = InequalityConstraint(Impact(link9,[0.0;.20;01.0]),Impact(link9,[0.0;-.20;01.0]))

# joint1 = InequalityConstraint(link1,cfr)
# joint2 = InequalityConstraint(link2,cfr)
# joint3 = InequalityConstraint(link3,cfr)
# joint4 = InequalityConstraint(link4,cfr)
# joint5 = InequalityConstraint(link5,cfr)
# joint6 = InequalityConstraint(link6,cfr)
# joint7 = InequalityConstraint(link7,cfr)
# joint8 = InequalityConstraint(link8,cfr)
# joint9 = InequalityConstraint(link9,cfr)
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
# link1.q[2] = Quaternion(AngleAxis(-rand()-0.2,rand(3)-ones(3)*0.5...))
link1.q[2] = Quaternion(SVector([0.9127362490430289;-0.19667611676218716;-0.01932481073253633;-0.3575718060311244]...))
# link1.q[2] = Quaternion(SVector([0.8264247451585373;0.007903070296695762;0.42943709893585524;-0.364065186645316]...))
# link1.q[2] = Quaternion(SVector([0.9011767166792705;0.037092558277805214;0.22503618549822949;-0.36859650385209014]...))
qtemp = link1.q[2]

for link in links
    # link.x[2] += [0.;0.05;0.05]
    link.x[2] += [0.0;0.03;0.07]
end

# simulate_ip!(mech,save=true,debug=false)
# MaximalCoordinateDynamics.visualize(mech,shapes)
