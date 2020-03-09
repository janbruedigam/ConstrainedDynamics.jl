using Rotations
using Plots: RGBA
using StaticArrays

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 1.0
x,z = .1,.1
b1 = Box(x,l1,z,l1*1000,color=RGBA(1.,1.,0.))
b2 = Box(l1,x,z,l1,color=RGBA(1.,0.,0.))
b3 = Box(l1*2,x,z,l1,color=RGBA(1.,0.,0.))
b4 = Box(x,l1,z,l1*1,color=RGBA(1.,0.,0.))

vert11 = [0.;l1/2;0.]
vert12 = -vert11

vert21 = [l1/2;0.;0.]
vert22 = -vert21


# Links
origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b2)
link3 = Body(b2)
link4 = Body(b4)

# Constraints
joint0to1 = EqualityConstraint(OriginConnection(origin,link1))
# joint0to1 = EqualityConstraint(Fixed(origin,link1,zeros(3),zeros(3)))
joint1to23 = EqualityConstraint(Fixed(link1,link2,vert12,vert22),Fixed(link1,link3,vert11,vert22))
joint1to4 = EqualityConstraint(Revolute(link1,link4,vert12,vert11,ex))



# links = [link1; link2; link3]
# constraints = [joint0to1;joint1to23]
links = [link1; link4]
constraints = [joint0to1;joint1to4]
shapes = [b1;b2;b3;b4]

mech = Mechanism(origin,links,constraints, shapes=shapes,g=0.)
# setPosition!(mech,link1,link2,p1=vert12,p2=vert22)
# setPosition!(mech,link1,link3,p1=vert11,p2=vert22)
# setForce!(mech,link3,τ=[0.1;0.;0.])
# setForce!(mech,link2,τ=[0;0;2.])
# setForce!(mech,link3,τ=[0;0;2.])
setPosition!(mech,link1,link4,p1=vert12,p2=vert11)
# setForce!(mech,link1,τ=[-0.05;0;0])
# setForce!(mech,link4,τ=[0.05;0;0])

setForce!([nothing;[0.05]],joint1to4,mech)

simulate!(mech,save=true)
visualize!(mech)
