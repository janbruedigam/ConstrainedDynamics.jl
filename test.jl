using ConstrainedDynamics

# Parameters
ex = [1.;0.;0.]

h = 1.
r = .05
b1 = Cylinder(r, h, h, color = RGBA(1., 0., 0.))

vert11 = [0.;0.;h / 2]
vert12 = -vert11

# Initial orientation
phi = pi / 4
q1 = Quaternion(RotX(phi))

# Links
N = 6

origin = Origin{Float64}()
links = [Body(b1) for i = 1:N]

# Constraints
cf = 0.2
ineqcs = [InequalityConstraint(Friction(links[i], [0;0;1.0],cf)) for i = 2:N-1]

jointb1 = EqualityConstraint(OriginConnection(origin, links[1]))
constraints = [jointb1;[EqualityConstraint(Revolute(links[i - 1], links[i], vert12, vert11, ex)) for i = 2:N]]

shapes = [b1]

mech = Mechanism(origin, links, constraints, ineqcs;tend = 20.,Δt = 0.01, shapes = shapes)
setPosition!(mech,origin,links[1],p2 = vert11, Δx=[0;0;N*1.], Δq = q1)
previd = links[1].id
for body in Iterators.drop(mech.bodies, 1)
    global previd
    setPosition!(mech, ConstrainedDynamics.getbody(mech, previd), body, p1 = vert12, p2 = vert11)
    previd = body.id
end

simulate!(mech,save = true,debug=false)
visualize!(mech)

# include("examples/n_pendulum.jl")
# ConstrainedDynamics.setentries!(mech)
# A=ConstrainedDynamics.formMatrix(mech)
# cond(A)