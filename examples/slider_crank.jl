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
setInit!(origin,link2,p3+p1,zeros(3),τ=[0;0;0.3])

# Constraints
joint0to12 = Constraint(CylindricalFree(origin,link1,zeros(3),p0,ey),Revolute(origin,link2,p3+p1,zeros(3),ez))
joint1to2 = Constraint(Cylindrical(link1,link2,p1,p3,ez))


links = [link1; link2]
constraints = [joint0to12;joint1to2]
shapes = [box,cyl]

mech = Mechanism(origin,links,constraints)
# link2.q[2] = Quaternion(RotZ(0.05))

simulate!(mech,save=true)
MaximalCoordinateDynamics.visualize(mech,shapes)
