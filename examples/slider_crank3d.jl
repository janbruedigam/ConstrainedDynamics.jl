using Rotations
using Plots: RGBA

!(@isdefined MaximalCoordinateDynamics) && include(joinpath("..", "src", "MaximalCoordinateDynamics.jl"))
using Main.MaximalCoordinateDynamics

# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 2.0
x,y = .1,.1
r = l1/4
h = r/10
box = Box(x,l1,y,l1,color=RGBA(1.,1.,0.))
cyl = Cylinder(r,h,2*r,color=RGBA(1.,0.,0.))
box2 = Box(h,2r,2r,2r,color=RGBA(1.,0.,1.))

p0 = [0;l1/2;0]
p1 = -p0
p2 = [0;r;0]
p3 = -p2

# Initial orientation
q1 = Quaternion(RotY(pi/2))

# Links
origin = Origin{Float64}()

link1 = Body(box)
setInit!(origin,link1,zeros(3),p0)

link2 = Body(cyl)
setInit!(origin,link2,p3+p1,zeros(3),τ=[0;0;1.])

link3 = Body(box2)
setInit!(origin,link3,p2,zeros(3))

# Constraints
joint0to23 = Constraint(Revolute(origin,link2,p3+p1,zeros(3),ez),Revolute(origin,link3,p2,zeros(3),ex))
joint1to23 = Constraint(Spherical(link1,link2,p1,p3),Spherical(link1,link3,p0,p3))



links = [link1;link2;link3]
constraints = [joint0to23;joint1to23]
shapes = [box,cyl,box2,indicator]

mech = Mechanism(origin,links,constraints,g=0.,tend=20.)
link2.q[2] = Quaternion(RotZ(0.01))
link3.q[2] = Quaternion(RotX(0.01))

simulate!(mech,save=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
