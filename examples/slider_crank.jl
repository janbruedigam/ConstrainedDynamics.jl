using StaticArrays
using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 2.0
x, y = .1, .1
r = l1 / 4
h = r / 10
box = Box(x, l1, y, l1, color = RGBA(1., 1., 0.))
cyl = Cylinder(r, h, 2 * r, color = RGBA(1., 0., 0.))

p0 = [0;l1 / 2;0]
p1 = -p0
p2 = [0;r;0]
p3 = -p2

# Initial orientation
q1 = Quaternion(RotY(pi / 2))

# Links
origin = Origin{Float64}()
link1 = Body(box)
link2 = Body(cyl)

# Constraints
joint0to12 = EqualityConstraint(CylindricalFree(origin, link1, zeros(3), p0, ey), Revolute(origin, link2, p3 + p1, zeros(3), ez))
joint1to2 = EqualityConstraint(Cylindrical(link1, link2, p1, p3, ez))


links = [link1; link2]
constraints = [joint0to12;joint1to2]
shapes = [box,cyl]


mech = Mechanism(origin, links, constraints, shapes = shapes)
setPosition!(mech,link1,x = -p0)
setPosition!(mech,link2,x = p3 + p1)

a1 = nothing # SVector{0,Float64}() # SVector{1,Float64}(0)
a2 = nothing # SVector{0,Float64}() # SVector{3,Float64}(0,0,0)
a3 = nothing # SVector{0,Float64}() # SVector{0,Float64}()
a4 = 0.3 # SVector{1,Float64}(0.3)

a1 = SVector{0,Float64}() # SVector{1,Float64}(0)
a2 = SVector{0,Float64}() # SVector{3,Float64}(0,0,0)
a3 = SVector{0,Float64}() # SVector{0,Float64}()
a4 = 0.3 # SVector{1,Float64}(0.3)

a = [[a1];[a2];[a3];[a4]]
# a = [a1;a2;a3;a4]

setForce!(a,joint0to12,mech)
# setForce!(mech,link2,τ=[0;0;0.3])

simulate!(mech,save = true)
visualize!(mech)
