using Rotations
using Plots: RGBA
using StaticArrays

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 0.5
width,depth = 0.5, 0.5
box1 = Box(width,depth,length1,1.,color=RGBA(1.,1.,0.))
b1 = Box(0.1,0.1,.1,.1,color=RGBA(1.,0.,0.))

# Links
origin = Origin{Float64}()

link1 = Body(box1)
link2 = Body(b1)
link3 = Body(b1)
link4 = Body(b1)
link5 = Body(b1)
link6 = Body(b1)
link7 = Body(b1)
link8 = Body(b1)
link9 = Body(b1)

# Constraints
joint1 = InequalityConstraint(Impact(link1,[0;0;1.0]))
joint2 = InequalityConstraint(Impact(link2,[0;0;1.0]))
joint3 = InequalityConstraint(Impact(link3,[0;0;1.0]))
joint4 = InequalityConstraint(Impact(link4,[0;0;1.0]))
joint5 = InequalityConstraint(Impact(link5,[0;0;1.0]))
joint6 = InequalityConstraint(Impact(link6,[0;0;1.0]))
joint7 = InequalityConstraint(Impact(link7,[0;0;1.0]))
joint8 = InequalityConstraint(Impact(link8,[0;0;1.0]))
joint9 = InequalityConstraint(Impact(link9,[0;0;1.0]))

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
shapes = [box1;b1]


mech = Mechanism(origin, links,constraints,ineqs,shapes=shapes)
setPosition!(mech,link1,x=[0.;-2;1.5])
setPosition!(mech,link1,link2,p1=[length1/2;length1/2;-length1/2])
setPosition!(mech,link1,link3,p1=[length1/2;-length1/2;-length1/2])
setPosition!(mech,link1,link4,p1=[-length1/2;length1/2;-length1/2])
setPosition!(mech,link1,link5,p1=[-length1/2;-length1/2;-length1/2])
setPosition!(mech,link1,link6,p1=[length1/2;length1/2;length1/2])
setPosition!(mech,link1,link7,p1=[length1/2;-length1/2;length1/2])
setPosition!(mech,link1,link8,p1=[-length1/2;length1/2;length1/2])
setPosition!(mech,link1,link9,p1=[-length1/2;-length1/2;length1/2])

setVelocity!(mech,link1,v=[0;3;7.],ω=(rand(3).-0.5)*100)
setVelocity!(mech,link1,link2,p1=[length1/2;length1/2;-length1/2])
setVelocity!(mech,link1,link3,p1=[length1/2;-length1/2;-length1/2])
setVelocity!(mech,link1,link4,p1=[-length1/2;length1/2;-length1/2])
setVelocity!(mech,link1,link5,p1=[-length1/2;-length1/2;-length1/2])
setVelocity!(mech,link1,link6,p1=[length1/2;length1/2;length1/2])
setVelocity!(mech,link1,link7,p1=[length1/2;-length1/2;length1/2])
setVelocity!(mech,link1,link8,p1=[-length1/2;length1/2;length1/2])
setVelocity!(mech,link1,link9,p1=[-length1/2;-length1/2;length1/2])


simulate!(mech,save=true)
visualize!(mech)
