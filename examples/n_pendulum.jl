using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05
b1 = Cylinder(r, h, h, color = RGBA(1., 0., 0.))

vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = UnitQuaternion(RotX(phi))

# Links
N = 20

origin = Origin{Float64}()
links = [Body(b1) for i = 1:N]

# Constraints
jointb1 = EqualityConstraint(Revolute(origin, links[1], ex; p2=vert11))
if N>1
    constraints = [jointb1;[EqualityConstraint(Revolute(links[i - 1], links[i], ex; p1=vert12, p2=vert11)) for i = 2:N]]
else
    constraints = [jointb1]
end

shapes = [b1]

mech = Mechanism(origin, links, constraints;Δt = 0.01, shapes = shapes)
setPosition!(origin,links[1],p2 = vert11,Δq = q1)
previd = links[1].id
for body in Iterators.drop(mech.bodies, 1)
    global previd
    setPosition!(ConstrainedDynamics.getbody(mech, previd), body, p1 = vert12, p2 = vert11)
    previd = body.id
end

storage = simulate!(mech, 10., record = true)
visualize(mech, storage, shapes)