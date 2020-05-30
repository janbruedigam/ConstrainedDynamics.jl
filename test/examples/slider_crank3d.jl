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
box2 = Box(h, 2r, 2r, 2r, color = RGBA(1., 0., 1.))

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
link3 = Body(box2)

# Constraints
joint0to23 = EqualityConstraint(Revolute(origin, link2, ez; p1=p3 + p1, p2=zeros(3)), Revolute(origin, link3, ex; p1=p2, p2=zeros(3)))
joint1to23 = EqualityConstraint(Spherical(link1, link2; p1=p1, p2=p3), Spherical(link1, link3; p1=p0, p2=p3))


links = [link1;link2;link3]
constraints = [joint0to23;joint1to23]
shapes = [box,cyl,box2]

mech = Mechanism(origin, links, constraints, g = 0., shapes = shapes)
setPosition!(link1,x = -p0)
setPosition!(link2,x = p3 + p1)
setPosition!(link3,x = -p2)
setVelocity!(link2,ω = [0;0;1.])
setVelocity!(link3,ω = [1.;0;0])
setForce!(link2,τ = [0;0;1.])

