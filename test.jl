using Rotations
using Plots: RGBA
using StaticArrays

!(@isdefined MaximalCoordinateDynamics) && include(joinpath(pwd(), "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 1.0
x,y = .1,.2
b1 = Box(x,y,l1,l1,color=RGBA(1.,1.,0.))

vert11 = [0.;0;l1/2]
vert12 = -vert11


# Links
origin = Origin{Float64}()
link1 = Body(b1)
link2 = Body(b1)

# Constraints
joint0to1 = EqualityConstraint(Revolute(origin,link1,zeros(3),vert11,ex,offset=Quaternion(RotY(pi/4))))



links = [link1]
constraints = [joint0to1]
shapes = [b1]

mech = Mechanism(origin,links,constraints, shapes=shapes,tend=0.01)
setPosition!(mech,origin,link1,p2=vert11,Δq=Quaternion(RotY(pi/4))*Quaternion(RotX(pi/4)))


# simulate!(mech,save=true)
# visualize!(mech)
