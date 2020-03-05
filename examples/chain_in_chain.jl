using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
ex = [1.;0.;0.]

l1 = 1.0
l2 = 1.0#sqrt(2)/2
x,y = .1,.1
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))
b2 = Box(x,y,l2,l2,color=RGBA(1.,1.,0.))
b3 = Box(x,y,l1,l1,color=RGBA(1.,0.,0.))
b4 = Box(x,y,l2,l2,color=RGBA(1.,0.,0.))

vert11 = [0.;0.;l1/2]
vert12 = -vert11

vert21 = [0.;0.;l2/2]
vert22 = -vert21

# Initial orientation
phi1 = pi/5
q1 = Quaternion(RotX(phi1))

# Links
origin = Origin{Float64}()

link1 = Body(b1)
link2 = Body(b2)
link3 = Body(b3)
link4 = Body(b4)

link5 = Body(b1)
link6 = Body(b2)
link7 = Body(b3)
link8 = Body(b4)

# Constraints
joint0to1 = EqualityConstraint(Revolute(origin,link1,zeros(3),vert11,ex))
joint1to23 = EqualityConstraint(Revolute(link1,link2,vert12,vert21,ex),Cylindrical(link1,link3,vert11,vert11,ex))
joint3to4 = EqualityConstraint(Revolute(link3,link4,vert12,vert21,ex))
joint2to4 = EqualityConstraint(Revolute(link2,link4,vert22,vert22,ex))

joint4to5 = EqualityConstraint(Revolute(link4,link5,zeros(3),vert11,ex))
joint5to67 = EqualityConstraint(Revolute(link5,link6,vert12,vert21,ex),Cylindrical(link5,link7,vert11,vert11,ex))
joint7to8 = EqualityConstraint(Revolute(link7,link8,vert12,vert21,ex))
joint6to8 = EqualityConstraint(Revolute(link6,link8,vert22,vert22,ex))


links = [link1; link2; link3; link4; link5; link6; link7; link8]
constraints = [joint0to1; joint1to23; joint3to4; joint2to4; joint4to5;joint5to67;joint7to8;joint6to8]
shapes = [b1,b2,b3,b4]

mech = Mechanism(origin,links, constraints,InequalityConstraint{Float64}[],shapes=shapes)

setPosition!(mech,origin,link1,p2=vert11,Δq=q1)
setPosition!(mech,link1,link2,p1=vert12,p2=vert21,Δq=inv(q1))
setPosition!(mech,link1,link3,p1=vert11,p2=vert11,Δq=inv(q1))
setPosition!(mech,link3,link4,p1=vert12,p2=vert21,Δq=q1)

setPosition!(mech,link4,link5,p2=vert11)
setPosition!(mech,link5,link6,p1=vert12,p2=vert21,Δq=inv(q1)*inv(q1))
setPosition!(mech,link5,link7,p1=vert11,p2=vert11,Δq=inv(q1)*inv(q1))
setPosition!(mech,link7,link8,p1=vert12,p2=vert21,Δq=q1*q1)


simulate!(mech,save=true)
visualize!(mech)
